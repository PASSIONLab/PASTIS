// Created by Saliya Ekanayake on 2019-09-03.

#include <algorithm>
#include <thread>
#include <utility>

#include <cuda_runtime_api.h>
#include <cuda.h>

#include "../../inc/util.hpp"
#include "../../inc/align/AdeptBSWAligner.hpp"

using std::string;	using std::vector;	using std::to_string;	using std::max;
using std::min;	using std::thread;

extern shared_ptr<pastis::ParallelOps> parops;

#define SANITIZE_PAIRS



namespace
pastis
{




void
AdeptBSWAligner::construct_seqs
(
    std::shared_ptr<DistFastaData> dfd
)
{
	uint64_t rseq_cnt = dfd->get_nrseqs();
	uint64_t cseq_cnt = dfd->get_ncseqs();

	rseqs_.resize(rseq_cnt);
	cseqs_.resize(cseq_cnt);

	int			cur_nbr	 = 0;
	uint64_t	cur_rseq = 0;
	uint64_t	cur_cseq = 0;
	for (auto &nbr : dfd->my_nbrs)
	{
		uint64_t nbr_seqs_count =
			nbr.nbr_seq_end_idx-nbr.nbr_seq_start_idx+1;
		FastaData	*curfd = dfd->fd;
		uint64_t	 sidx  = nbr.nbr_seq_start_idx;
		if (nbr.nbr_rank != parops->g_rank) // not myself
		{
			curfd = dfd->recv_fds_[cur_nbr];
			sidx  = 0;
			++cur_nbr;
		}		

		vector<string> &cur_seqs = (nbr.rc_flag==1 ? rseqs_ : cseqs_);
		uint64_t seq_beg = (nbr.rc_flag==1 ? cur_rseq : cur_cseq);

		#pragma omp parallel
		{
			ushort len;
			uint64_t start_offset, end_offset_inclusive;

			#pragma omp for
			for (uint64_t i = 0; i < nbr_seqs_count; ++i)
			{
				char *buff =
					curfd->get_sequence(sidx+i, len, start_offset,
										end_offset_inclusive);
				cur_seqs[seq_beg+i] = std::move(string(buff+start_offset, len));
			}
		}

		if (nbr.rc_flag == 1)
			cur_rseq += nbr_seqs_count;
		else
			cur_cseq += nbr_seqs_count;		
	}


	// diag cell row and col seqs are the same
	if (parops->grid->GetRankInProcRow() == parops->grid->GetRankInProcCol())
	{
		assert(cseqs_.empty());
		cseqs_.assign(rseqs_.begin(), rseqs_.end());
	}
}



void
AdeptBSWAligner::construct_seqs_bl
(
    std::shared_ptr<DistFastaData> dfd
)
{
	string s_tmp;
	parops->logger->log("constructing library-specific seqs.");
	
	uint64_t rseq_cnt = dfd->get_nrseqs_bl();
	uint64_t cseq_cnt = dfd->get_ncseqs_bl();

	parops->logger->log("#rseqs " + to_string(rseq_cnt) +
						" #cseqs " + to_string(cseq_cnt));


	rseqs_.resize(rseq_cnt);
	cseqs_.resize(cseq_cnt);

	int			cur_nbr	 = 0;
	uint64_t	cur_rseq = 0;
	uint64_t	cur_cseq = 0;
	for (auto &nbr : dfd->bl_nbrs_)
	{
		uint64_t nbr_seqs_count =
			nbr.nbr_seq_end_idx-nbr.nbr_seq_start_idx+1;
		FastaData	*curfd = dfd->fd;
		uint64_t	 sidx  = nbr.nbr_seq_start_idx;
		if (nbr.nbr_rank != parops->g_rank) // not myself
		{
			curfd = dfd->bl_rbufs_[cur_nbr].second;
			sidx  = 0;
			++cur_nbr;
		}

		// s_tmp = "nbr rank " + to_string(nbr.nbr_rank) +
		// 	" count " + to_string(nbr_seqs_count) +
		// 	" row? " + to_string(nbr.rc_flag==1) +
		// 	" cur_rseq " + to_string(cur_rseq) +
		// 	" cur_cseq " + to_string(cur_cseq) +
		// 	" fd addr " + ptr_to_str(curfd) +
		// 	" fd l_start " + to_string(curfd->get_l_start()) +
		// 	" fd l_end " + to_string(curfd->get_l_end());
		// parops->logger->log(s_tmp);

		vector<string> &cur_seqs = (nbr.rc_flag==1 ? rseqs_ : cseqs_);
		uint64_t seq_beg = (nbr.rc_flag==1 ? cur_rseq : cur_cseq);

		#pragma omp parallel
		{
			ushort len;
			uint64_t start_offset, end_offset_inclusive;

			#pragma omp for
			for (uint64_t i = 0; i < nbr_seqs_count; ++i)
			{
				char *buff =
					curfd->get_sequence(sidx+i, len, start_offset,
										end_offset_inclusive);

				string seq(buff+start_offset, len);
				cur_seqs[seq_beg+i] = std::move(seq);
			}
		}

		if (nbr.rc_flag == 1)
			cur_rseq += nbr_seqs_count;
		else
			cur_cseq += nbr_seqs_count;		
	}

	parops->logger->log("completed constructing sequences.");
}




void
AdeptBSWAligner::aln_batch
(
    std::tuple<uint64_t, uint64_t, CommonKmerLight *>	*mattuples,
	uint64_t											 beg,
	uint64_t											 end,
	uint64_t 											 bl_roffset,
	uint64_t 											 bl_coffset,
	const params_t										&params	
)
{
	// parops->info("Adept BSW Aligner aln_batch");
	parops->logger->log("Adept BSW aln_batch beginning.");
	
	parops->tp->start_timer("align|pre");

	uint64_t		npairs		= end-beg;
	vector<string>	seqs_q(npairs);	// queries - shorter seqs
	vector<string>	seqs_r(npairs);	// refs - longer seqs
	uint64_t		max_rlen = 0;
	uint64_t		max_qlen = 0;
	// ADEPT needs uppercase seqs
	auto f_upper = []
		(char c)
		{
			return static_cast<char>(std::toupper(c));
		};

	int numThreads = 1;
	#ifdef THREADED
	#pragma omp parallel
    {
      	numThreads = omp_get_num_threads();
    }
	#endif

	parops->logger->log("Preparing sequences for ADEPT");

	uint64_t mem_cur, mem_tot = 0, interval = 1e6;
	double aln_lens = 0;

	parops->tp->mat_stats["aln_pairs"] += npairs;
	
	#pragma omp parallel for reduction(max: max_rlen, max_qlen) reduction(+: aln_lens)
	for (uint64_t i = beg; i < end; ++i)
	{
		uint64_t lr = std::get<0>(mattuples[i]);
		uint64_t lc = std::get<1>(mattuples[i]);
		assert(lr+bl_roffset < rseqs_.size() &&
			   lc+bl_coffset < cseqs_.size());
		string &rseq = rseqs_[lr+bl_roffset];
		string &cseq = cseqs_[lc+bl_coffset];
		if (rseq.size() < cseq.size())
		{
			seqs_q[i-beg] = rseq;
			seqs_r[i-beg] = cseq;
		}
		else
		{
			seqs_q[i-beg] = cseq;
			seqs_r[i-beg] = rseq;
		}

		std::transform(seqs_q[i-beg].begin(), seqs_q[i-beg].end(),
					   seqs_q[i-beg].begin(), f_upper);
		std::transform(seqs_r[i-beg].begin(), seqs_r[i-beg].end(),
					   seqs_r[i-beg].begin(), f_upper);
		max_rlen = max(max_rlen, seqs_r[i-beg].size());
		max_qlen = max(max_qlen, seqs_q[i-beg].size());

		aln_lens += seqs_r[i-beg].size() * seqs_q[i-beg].size();
	}

	parops->tp->mat_stats["aln_pair_lens"] += aln_lens;

	parops->tp->stop_timer("align|pre");

	parops->tp->start_timer("align|multi_gpu");

	parops->logger->log("Adept BSW aln_batch multi_gpu function call.");

	auto t_begin = std::chrono::system_clock::now();

	auto all_results =
		std::move(ADEPT::multi_gpu(seqs_r, seqs_q,
						 ADEPT::options::ALG_TYPE::SW,
						 ADEPT::options::SEQ_TYPE::AA,
						 ADEPT::options::CIGAR::NO,
						 ADEPT::options::SCORING::ALNS_AND_SCORE,
						 max_rlen, max_qlen,
						 score_mat_, gaps_, g_batch_sz_));
	
	parops->tp->stop_timer("align|multi_gpu");

	parops->logger->log("Adept BSW aln_batch alignments complete.");

	omp_set_num_threads(numThreads);

	parops->tp->start_timer("align|post");

	uint64_t cur_beg = 0;
	for (int gpu_idx = 0; gpu_idx < all_results.gpus; ++gpu_idx)
	{
		uint64_t cur_cnt = all_results.per_gpu;
		if (gpu_idx == all_results.gpus-1)
            cur_cnt += all_results.left_over;

		auto &res = all_results.results[gpu_idx];
		
		#pragma omp for
		for (uint64_t i = 0; i < cur_cnt; ++i)
		{
			uint64_t cur = cur_beg+i;
			int len_r = seqs_r[cur].size();
			int len_q = seqs_q[cur].size();
			int alen_r = res.ref_end[i]-res.ref_begin[i];
			int alen_q = res.query_end[i]-res.query_begin[i];

			double cov_r = (double)(alen_r) / len_r;
			double cov_q = (double)(alen_q) / len_q;

			if (max(cov_r, cov_q) >=
				params.aln_cov_thr) // coverage constraint
			{
				CommonKmerLight *ckl = std::get<2>(mattuples[beg+cur]);
				ckl->score_aln = (float)(res.top_scores[i]) /
					(float)min(len_r, len_q);
			}
		}
		cur_beg += cur_cnt;
	}

	parops->logger->log("Adept BSW aln_batch post-processing complete.");
	parops->tp->stop_timer("align|post");
}



void
AdeptBSWAligner::aln_batch_ovlp
(
    std::tuple<uint64_t, uint64_t, CommonKmerLight *>	*mattuples,
	uint64_t											 beg,
	uint64_t											 end,
	uint64_t 											 bl_roffset,
	uint64_t 											 bl_coffset,
	MultData<MatrixEntry, CommonKmerLight>				&md,
	const params_t										&params	
)
{
	// parops->info("Adept BSW Aligner aln_batch");
	parops->logger->log("Adept BSW aln_batch beginning.");
	
	parops->tp->start_timer("align|pre");

	uint64_t		npairs		= end-beg;
	vector<string>	seqs_q(npairs);	// queries - shorter seqs
	vector<string>	seqs_r(npairs);	// refs - longer seqs
	uint64_t		max_rlen = 0;
	uint64_t		max_qlen = 0;
	// ADEPT needs uppercase seqs
	auto f_upper = []
		(char c)
		{
			return static_cast<char>(std::toupper(c));
		};

	int numThreads = 1;
	#ifdef THREADED
	#pragma omp parallel
    {
      	numThreads = omp_get_num_threads();
    }
	#endif

	parops->logger->log("Preparing sequences for ADEPT");

	if (mattuples == nullptr && end != 0)
		parops->logger->log("mattuples empty but end index != 0!");

	uint64_t mem_cur, mem_tot = 0, interval = 1e6;
	double aln_lens = 0;

	parops->tp->mat_stats["aln_pairs"] += npairs;
	
	#pragma omp parallel for reduction(max: max_rlen, max_qlen) reduction(+: aln_lens)
	for (uint64_t i = beg; i < end; ++i)
	{
		uint64_t lr = std::get<0>(mattuples[i]);
		uint64_t lc = std::get<1>(mattuples[i]);
		assert(lr+bl_roffset < rseqs_.size() &&
			   lc+bl_coffset < cseqs_.size());
		string &rseq = rseqs_[lr+bl_roffset];
		string &cseq = cseqs_[lc+bl_coffset];
		if (rseq.size() < cseq.size())
		{
			seqs_q[i-beg] = rseq;
			seqs_r[i-beg] = cseq;
		}
		else
		{
			seqs_q[i-beg] = cseq;
			seqs_r[i-beg] = rseq;
		}

		std::transform(seqs_q[i-beg].begin(), seqs_q[i-beg].end(),
					   seqs_q[i-beg].begin(), f_upper);
		std::transform(seqs_r[i-beg].begin(), seqs_r[i-beg].end(),
					   seqs_r[i-beg].begin(), f_upper);
		max_rlen = max(max_rlen, seqs_r[i-beg].size());
		max_qlen = max(max_qlen, seqs_q[i-beg].size());

		aln_lens += seqs_r[i-beg].size() * seqs_q[i-beg].size();
		
		// if (i % interval == 0)
		// {
		// 	// mem_cur = 2 * interval * sizeof(string); // base strings
		// 	// for (uint64_t j = i-interval; j < i; ++j)
		// 	// 	mem_cur += seqs_q[j-beg].capacity() + seqs_r[j-beg].capacity();
		// 	// mem_tot += mem_cur;
		// 	parops->logger->log("Seq " + to_string(i));
		// 	// parops->logger->log("String pairs mem cur " + gb_str(mem_cur) +
		// 	// 					" total " + gb_str(mem_tot));			
		// }
	}

	parops->tp->mat_stats["aln_pair_lens"] += aln_lens;

	////////////////////////////////////////////////////////////////////////////
	// mem_tot = 2 * npairs * sizeof(string);
	// for (uint64_t i = beg; i < end; ++i)
	// 	mem_tot += seqs_q[i-beg].capacity() + seqs_r[i-beg].capacity();
	// parops->logger->log("String pairs " + to_string(npairs) +
	// 					" mem " + gb_str(mem_tot));
	////////////////////////////////////////////////////////////////////////////
	
	parops->tp->stop_timer("align|pre");

	parops->logger->log("Adept BSW aln_batch multi_gpu function call.");

	auto t_begin_02 = std::chrono::system_clock::now();
	parops->logger->log("Pre-align complete, will perform mult & alignment");

	int ngpus;
	cudaGetDeviceCount(&ngpus);

	parops->tp->start_timer("align|ovlp");
	
	thread thd_mult, thd_aln;
	bool is_mult_active = false;
	if ((md.br_ * md.bc_ > 1)	// first block computed earlier
		&& !(md.bri_ == md.br_-1 && md.bci_ == md.bc_-1) // not last block
		)
	{
		is_mult_active = true;
		int bri, bci;

		
		if (params.lb == params_t::LoadBal::LB_IDX)
		{
			parops->info("Before calling computeNext "  +
						 to_string(md.bri_) + " " + to_string(md.bci_));
			parops->logger->log("Creating multiplication thread");
			thd_mult = thread{&MultData<MatrixEntry,
							  CommonKmerLight>::computeNext,
							  &md, numThreads-6};
		}
		else
		{
			parops->info("Before calling computeNext_trg "  +
						 to_string(md.bri_) + " " + to_string(md.bci_));
			parops->logger->log("Creating multiplication thread");
			thd_mult = thread{&MultData<MatrixEntry,
							  CommonKmerLight>::computeNext_trg,
							  &md, numThreads-6};
		}
	}

	auto t_begin = std::chrono::system_clock::now();

	parops->logger->log("Calling multi_gpu");

	parops->tp->start_timer("align|multi_gpu");

	ADEPT::all_alns all_results(0);
	if (mattuples != nullptr)
		all_results = std::move
			(ADEPT::multi_gpu(seqs_r, seqs_q,ADEPT::options::ALG_TYPE::SW,
							  ADEPT::options::SEQ_TYPE::AA,
							  ADEPT::options::CIGAR::NO,
							  ADEPT::options::SCORING::ALNS_AND_SCORE,
							  max_rlen, max_qlen,
							  score_mat_, gaps_, g_batch_sz_));

	parops->tp->stop_timer("align|multi_gpu");
	parops->tp->start_timer("align|join-wait");

	if (is_mult_active)
	{
		parops->logger->log("Before joining with mult thread");
		thd_mult.join();
		parops->logger->log("joined with multiplication thread " +
					 to_string(md.bri_) + " " + to_string(md.bci_));
	}

	parops->tp->stop_timer("align|join-wait");

	double t_mult_align = 
			ms_t(std::chrono::system_clock::now()-t_begin_02).count();
	parops->info("Mult+align took " + to_string(t_mult_align/1e3) + " seconds");

	parops->logger->log("Adept BSW aln_batch alignments complete.");

	parops->tp->stop_timer("align|ovlp");

	omp_set_num_threads(numThreads);

	parops->tp->start_timer("align|post");

	if (mattuples != nullptr)
	{
		uint64_t cur_beg = 0;
		for (int gpu_idx = 0; gpu_idx < all_results.gpus; ++gpu_idx)
		{
			uint64_t cur_cnt = all_results.per_gpu;
			if (gpu_idx == all_results.gpus-1)
				cur_cnt += all_results.left_over;

			auto &res = all_results.results[gpu_idx];

			#pragma omp for
			for (uint64_t i = 0; i < cur_cnt; ++i)
			{
				uint64_t cur = cur_beg+i;
				int len_r = seqs_r[cur].size();
				int len_q = seqs_q[cur].size();
				int alen_r = res.ref_end[i]-res.ref_begin[i];
				int alen_q = res.query_end[i]-res.query_begin[i];

				double cov_r = (double)(alen_r) / len_r;
				double cov_q = (double)(alen_q) / len_q;

				if (max(cov_r, cov_q) >=
					params.aln_cov_thr) // coverage constraint
				{
					CommonKmerLight *ckl = std::get<2>(mattuples[beg+cur]);
					ckl->score_aln = (float)(res.top_scores[i]) /
						(float)min(len_r, len_q);
				}
			}
			cur_beg += cur_cnt;
		}
	}
	
	parops->logger->log("Adept BSW aln_batch post-processing complete.");

	parops->tp->stop_timer("align|post");
}


	
}

