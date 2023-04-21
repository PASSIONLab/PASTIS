// Created by Saliya Ekanayake on 2019-09-03.

#include <string>
#include <utility>

#include "seqan/align.h"
#include "seqan/align_parallel.h"
#include "seqan/basic.h"
#include "seqan/seeds.h"
#include "seqan/sequence.h"
#include "seqan/stream.h"

#include "../../inc/util.hpp"
#include "../../inc/align/SeqanXdropAligner.hpp"

using std::string;	using std::vector;	using std::to_string;

extern shared_ptr<pastis::ParallelOps> parops;

typedef seqan::Seed<seqan::Simple> TSeed;




namespace
pastis
{




void
SeqanXdropAligner::construct_seqs
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

		vector<seqan::Peptide> &cur_seqs = (nbr.rc_flag==1 ? rseqs_ : cseqs_);
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
				cur_seqs[seq_beg+i] =
					std::move(seqan::Peptide(buff+start_offset, len));
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

	#if PASTIS_DBG_LVL > 0
	uint64_t tot_sz = 0;
	for (auto &el : rseqs_)
		tot_sz += seqan::capacity(el) * sizeof(char);
	for (auto &el : cseqs_)
		tot_sz += seqan::capacity(el) * sizeof(char);
	parops->logger->log("approximate memory allocated for seqan seqs " +
						gb_str(tot_sz));
	parops->bytes_alloc += tot_sz;
	parops->logger->log("approximate memory in usage " +
						gb_str(parops->bytes_alloc));
	#endif
}



void
SeqanXdropAligner::construct_seqs_bl
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

		vector<seqan::Peptide> &cur_seqs = (nbr.rc_flag==1 ? rseqs_ : cseqs_);
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

				// @OGUZ-WARNING Seqan fails below
				// cur_seqs[seq_beg+i] =
				// 	std::move(seqan::Peptide(buff+start_offset, len));				
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
SeqanXdropAligner::aln_batch
(
    std::tuple<uint64_t, uint64_t, CommonKmerLoc *>		*mattuples,
	uint64_t											 beg,
	uint64_t											 end,
	uint64_t 											 bl_roffset,
	uint64_t 											 bl_coffset,
	const params_t										&params	
)
{
	parops->tp->start_timer("align|pre");
	uint64_t npairs = end-beg;
		
	// form seqan pairs
	seqan::StringSet<seqan::Gaps<seqan::Peptide>> seqsr;
	seqan::StringSet<seqan::Gaps<seqan::Peptide>> seqsc;
	resize(seqsr, npairs, seqan::Exact{});
	resize(seqsc, npairs, seqan::Exact{});

	parops->tp->mat_stats["aln_pairs"] += npairs;

	double aln_lens = 0;

	#pragma omp parallel for reduction(+: aln_lens)
	for (uint64_t i = beg; i < end; ++i)
	{
		uint64_t lr = std::get<0>(mattuples[i]);
		uint64_t lc = std::get<1>(mattuples[i]);
		assert(lr+bl_roffset < rseqs_.size() &&
			   lc_bl_coffset < cseqs_.size());
		
		seqsr[i-beg] = std::move(seqan::Gaps<seqan::Peptide>
								 (rseqs_[lr+bl_roffset]));
		seqsc[i-beg] = std::move(seqan::Gaps<seqan::Peptide>
								 (cseqs_[lc+bl_coffset]));

		aln_lens += seqan::length(seqan::source(seqsc[i-beg])) +
							   seqan::length(seqan::source(seqsr[i-beg]));
	}

	parops->tp->mat_stats["aln_pair_lens"] += aln_lens;
	parops->tp->stop_timer("align|pre");

	#if PASTIS_DBG_LVL > 0
	uint64_t tot_sz = 0;
	for (auto &el : seqsr)
		tot_sz += seqan::capacity(el) * sizeof(char);
	for (auto &el : seqsc)
		tot_sz += seqan::capacity(el) * sizeof(char);
	parops->logger->log("approximate memory allocated for alignment seqs " +
						gb_str(tot_sz));
	// parops->bytes_alloc += tot_sz;
	// parops->logger->log("approximate memory in usage " +
	// 					gb_str(parops->bytes_alloc));
	#endif

	// parallel alignment delegated to seqan
	int nthds = 1;
	#ifdef THREADED
	#pragma omp parallel
    {
      	nthds = omp_get_num_threads();
    }
	#endif
	
	seqan::ExecutionPolicy<seqan::Parallel, seqan::Vectorial> exec_policy;
	setNumThreads(exec_policy, nthds);

	std::vector<seqan::AlignmentStats> stats(npairs);
	for (int cnt = 0; cnt < seed_cnt_; ++cnt)
	{
		parops->tp->start_timer("align|pre");
		
		seqan::StringSet<seqan::Gaps<seqan::Peptide>> seqsc_ex;
		seqan::StringSet<seqan::Gaps<seqan::Peptide>> seqsr_ex;
		resize(seqsc_ex, npairs, seqan::Exact{});
		resize(seqsr_ex, npairs, seqan::Exact{});

		// extend the current seed and form a new gaps object
		#pragma omp parallel for
		for (uint64_t i = 0; i < npairs; ++i)
		{
			CommonKmerLoc *ckl = std::get<2>(mattuples[beg+i]);
			ushort l_row_seed_start_offset =
				(cnt == 0) ? ckl->first.first : ckl->second.first;
			ushort l_col_seed_start_offset =
				(cnt == 0) ? ckl->first.second : ckl->second.second;

			TSeed seed(l_col_seed_start_offset, l_row_seed_start_offset,
					   seed_len_);
			extendSeed(seed, seqan::source(seqsc[i]), seqan::source(seqsr[i]),
					   seqan::EXTEND_BOTH, blosum62_, xdrop_,
					   seqan::GappedXDrop());
			assignSource(seqsc_ex[i],
						 infix(seqan::source(seqsc[i]),
							   beginPositionH(seed), endPositionH(seed)));
			assignSource(seqsr_ex[i],
						 infix(seqan::source(seqsr[i]),
							   beginPositionV(seed), endPositionV(seed)));
		}

		parops->tp->stop_timer("align|pre");
		parops->tp->start_timer("align|cpu");

		// alignment
		globalAlignment(exec_policy, seqsc_ex, seqsr_ex, blosum62_);

		parops->tp->stop_timer("align|cpu");

		parops->tp->start_timer("align|post");

		// stats
		if (cnt == 0)
		{
			#pragma omp parallel for
			for (uint64_t i = 0; i < npairs; ++i)
			{
				computeAlignmentStats
					(stats[i], seqsc_ex[i], seqsr_ex[i], blosum62_);
			}
		}
		else
		{
			#pragma omp parallel for
			for (uint64_t i = 0; i < npairs; ++i)
			{
				seqan::AlignmentStats tmp;
				computeAlignmentStats(tmp, seqsc_ex[i], seqsr_ex[i],
									  blosum62_);
				if (tmp.alignmentIdentity > stats[i].alignmentIdentity)
					stats[i] = std::move(tmp);
			}
		}

		parops->tp->stop_timer("align|post");
	}
	

	parops->tp->start_timer("align|post");

	// stats
	#pragma omp parallel
	{		
		#pragma omp for
		for (uint64_t i = 0; i < npairs; ++i)
		{
			seqan::AlignmentStats &tmp = stats[i];
			double alen_minus_gapopens =
				tmp.alignmentLength - tmp.numGapOpens;			
			int len_seqc = seqan::length(seqan::source(seqsc[i]));
 			int len_seqr = seqan::length(seqan::source(seqsr[i]));

			// only keep alignments that meet coverage and ani criteria
			if (std::max((alen_minus_gapopens/len_seqc),
						 (alen_minus_gapopens/len_seqr)) >= params.aln_cov_thr
				&&
				tmp.alignmentIdentity >= params.aln_ani_thr)
			{
				CommonKmerLoc *ckl = std::get<2>(mattuples[beg+i]);
				ckl->score_aln =
					static_cast<float>(tmp.alignmentIdentity)/100.0f;
			}
		}
	}

	// #if PASTIS_DBG_LVL > 0
	// parops->bytes_alloc -= tot_sz;
	// #endif

	parops->tp->stop_timer("align|post");
}
	
}
