// Created by Saliya Ekanayake on 2019-09-03.

#include <string>
#include <utility>

#include "seqan/align.h"
#include "seqan/align_parallel.h"
#include "seqan/basic.h"
#include "seqan/stream.h"

#include "../../inc/util.hpp"
#include "../../inc/align/SeqanFullAligner.hpp"

using std::string;	using std::vector;	using std::to_string;

extern shared_ptr<pastis::ParallelOps> parops;




namespace
pastis
{




void
SeqanFullAligner::construct_seqs
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
}



void
SeqanFullAligner::construct_seqs_bl
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
SeqanFullAligner::aln_batch
(
    std::tuple<uint64_t, uint64_t, CommonKmerLight *>	*mattuples,
	uint64_t											 beg,
	uint64_t											 end,
	uint64_t 											 bl_roffset,
	uint64_t 											 bl_coffset,
	const params_t										&params	
)
{
	parops->tp->start_timer("sim:align_pre");
	
	// form seqan pairs
	seqan::StringSet<seqan::Gaps<seqan::Peptide>> seqsr;
	seqan::StringSet<seqan::Gaps<seqan::Peptide>> seqsc;
	resize(seqsr, end-beg, seqan::Exact{});
	resize(seqsc, end-beg, seqan::Exact{});

	#pragma omp for
	for (uint64_t i = beg; i < end; ++i)
	{
		uint64_t lr = std::get<0>(mattuples[i]);
		uint64_t lc = std::get<1>(mattuples[i]);
		assert(lr+bl_roffset < rseqs_.size() &&
			   lc+bl_coffset < cseqs_.size());
		
		seqsr[i-beg] = std::move(seqan::Gaps<seqan::Peptide>
								 (rseqs_[lr+bl_roffset]));
		seqsc[i-beg] = std::move(seqan::Gaps<seqan::Peptide>
								 (cseqs_[lc+bl_coffset]));
	}

	parops->tp->stop_timer("sim:align_pre");


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

	parops->tp->start_timer("sim:align");

	localAlignment(exec_policy, seqsc, seqsr, blosum62_);

	parops->tp->stop_timer("sim:align");


	parops->tp->start_timer("sim:align_post");

	// stats
	#pragma omp parallel
	{
		seqan::AlignmentStats stats;
		
		#pragma omp for
		for (uint64_t i = 0; i < end-beg; ++i)
		{
			computeAlignmentStats(stats, seqsc[i], seqsr[i], blosum62_);
			double alen_minus_gapopens =
				stats.alignmentLength - stats.numGapOpens;
			int len_seqc = seqan::length(seqan::source(seqsc[i]));
 			int len_seqr = seqan::length(seqan::source(seqsr[i]));

			// only keep alignments that meet coverage and ani criteria
			if (std::max((alen_minus_gapopens/len_seqc),
						 (alen_minus_gapopens/len_seqr)) >= params.aln_cov_thr
				&&
				stats.alignmentIdentity >= params.aln_ani_thr)
			{
				CommonKmerLight *ckl = std::get<2>(mattuples[beg+i]);
				ckl->score_aln =
					static_cast<float>(stats.alignmentIdentity)/100.0f;
			}
		}
	}

	parops->tp->stop_timer("sim:align_post");
}
	
}
