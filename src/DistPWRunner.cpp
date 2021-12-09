//
// Created by Saliya Ekanayake on 2019-02-19.
// Modified by Aydin Buluc on 2019-12-29 
//

#include <iostream>
#include <string>
#include <tuple>

#include "../inc/align/PWAlign.hpp"
#include "../inc/DistPWRunner.hpp"
#include "../inc/util.hpp"

using std::cout;	using std::endl;	using std::tuple;	using std::string;
using std::to_string;

extern shared_ptr<pastis::ParallelOps> parops;




namespace
pastis
{

template <typename NT>
void
pw_aln_batch
(
    std::shared_ptr<DistFastaData>		 dfd,
	combblas::SpParMat <uint64_t, NT, combblas::SpDCCols <uint64_t, NT>> &C,
	PWAlign 							&pwa,
	params_t							&params,
	int									 bri,
	int									 bci
)
{
	string		s_tmp;
	uint64_t	nrows = C.getnrow();
	uint64_t	ncols = C.getncol();
	uint64_t	nnzs  = C.getnnz();

	parops->info("Starting batch pairwise alignments.\n  "
				 "#elems prior to pre-pruning: " + to_string(nnzs));
	
	parops->tp->start_timer("sim:prune");
	
	// pre-prune the overlap matrix with symmetricity and kmer threshold
	// constraints
	// @OGUZ-TODO can do these both at once (kept for now for getting
	// eliminiation stats)
	if (bri == -1)				// not blocked
	{
		auto loc_up_tri =
			[&nrows, &ncols, &C] (const tuple<uint64_t, uint64_t, NT> &t)
			{
				uint64_t lr, lc;
				uint64_t gr = std::get<0>(t), gc = std::get<1>(t);
				C.Owner(nrows, ncols, gr, gc, lr, lc);
				return ((lc < lr) ||
						((lc == lr) && (gc <= gr)));			
			};
			C.PruneI(loc_up_tri);
	}
	
	uint64_t nnzs_sym = C.getnnz();
	parops->info("  #elems (sym. constraint): " +
				 to_string(nnzs_sym) + " (#elim " +
				 to_string(nnzs-nnzs_sym) + " " +
				 to_string((static_cast<double>(nnzs-nnzs_sym)/nnzs)*100.0)
				 + "%)");

	auto kmer_thr =
		[&params] (const NT &v)
		{
			return ((v.count <= params.ckthr) ||
					(v.score <= params.mosthr));
		};
	C.Prune(kmer_thr);

	uint64_t nnzs_thr = C.getnnz();

	// volatile int i = 0;
	// char hostname[256];
	// gethostname(hostname, sizeof(hostname));
	// printf("PID %d on %s ready for attach\n", getpid(), hostname);
	// fflush(stdout);
	// while (0 == i)
	// 	sleep(5);
	
	parops->info("  #elems (kmer thresholds): " +
				 to_string(nnzs_thr) + " (#elim " +
				 to_string(nnzs_sym-nnzs_thr) + " " +
				 to_string((static_cast<double>(nnzs_sym-nnzs_thr)/
							nnzs_sym)*100.0)
				 + "%)");

	
	// ADEPT aligner: both seqs cannot be greater than 1024
	if (params.pw_aln == params_t::PwAln::ALN_ADEPT_GPUBSW)
	{
		auto seq_len_filter =
			[&] (const tuple<uint64_t, uint64_t, NT> &t)
			{
				uint64_t lr, lc;
				uint64_t gr = std::get<0>(t), gc = std::get<1>(t);
				C.Owner(nrows, ncols, gr, gc, lr, lc);
				uint64_t bl_roffset = bri == -1 ? 0 : dfd->bl_rseq_offset(bri);
				uint64_t bl_coffset = bci == -1 ? 0 : dfd->bl_cseq_offset(bci);
				return
					(pwa.rseq_len(lr+bl_roffset) >= 1024) &&
					(pwa.cseq_len(lc+bl_coffset) >= 1024);
			};
			C.PruneI(seq_len_filter);

			uint64_t nnzs_seq_len = C.getnnz();
			parops->info("  #elems (seq len filter): " +
					 to_string(nnzs_seq_len) + " (#elim " +
					 to_string(nnzs_thr-nnzs_seq_len) + " " +
					 to_string((static_cast<double>(nnzs_thr-nnzs_seq_len)/
								nnzs_thr)*100.0)
					 + "%)");
	}

	parops->tp->stop_timer("sim:prune");

	parops->tp->start_timer("sim:align_all");

	// batch alignment
	// @OGUZ-TODO make this parallel
	uint64_t	l_nnz		= C.seqptr()->getnnz();	
	parops->logger->log("l_nnz before alignment " + to_string(l_nnz));	
	if (l_nnz != 0)
	{
		uint64_t	batch_sz	= params.aln_batch_sz;
		uint64_t	batch_cnt	= l_nnz/batch_sz + 1;
		uint64_t	batch_idx	= 0;
		
		// get tuples (need local indices)
		tuple<uint64_t, uint64_t, NT *> *mattuples =
			new tuple<uint64_t, uint64_t, NT *>[l_nnz];
		auto dcsc = C.seqptr()->GetDCSC();
		uint64_t k = 0;
		for (uint64_t i = 0; i < dcsc->nzc; ++i)
		{
			for (uint64_t j = dcsc->cp[i]; j < dcsc->cp[i+1]; ++j)
			{
				std::get<0>(mattuples[k]) = dcsc->ir[j];
				std::get<1>(mattuples[k]) = dcsc->jc[i];
				std::get<2>(mattuples[k]) = &(dcsc->numx[j]);
				++k;
			}
		}

		#if PASTIS_DBG_LVL > 0
		uint64_t tot_sz = (sizeof(uint64_t)*2 + sizeof(NT *)) * l_nnz;
		parops->bytes_alloc += tot_sz;
		parops->logger->log("approximate memory in usage " +
							gb_str(parops->bytes_alloc));
		#endif
		
		while (batch_idx < batch_cnt)
		{
			uint64_t beg = batch_idx * batch_sz;
			uint64_t end = ((batch_idx+1)*batch_sz > l_nnz) ?
				l_nnz : ((batch_idx+1)*batch_sz);
			parops->logger->log("processing batch " +
								to_string(batch_idx) + " of " +
								to_string(batch_cnt) +
								" size " + to_string(end-beg));
			pwa.aln_batch(mattuples, beg, end,
						  (bri == -1 ? 0 : dfd->bl_rseq_offset(bri)),
						  (bci == -1 ? 0 : dfd->bl_cseq_offset(bci)),
						  params);
			++batch_idx;
		}

		delete[] mattuples;

		#if PASTIS_DBG_LVL > 0
		parops->bytes_alloc -= tot_sz;
		#endif		
	}

	parops->tp->stop_timer("sim:align_all");

	parops->tp->start_timer("sim:prune");

	auto aln_prune = [] (NT &el)
		{
			return el.score_aln == 0.0f;
		};
	C.Prune(aln_prune);

	parops->tp->stop_timer("sim:prune");

	uint64_t nnzs_aln = C.getnnz();
	parops->info("  #elems (alignment thresholds): " +
				 to_string(nnzs_aln) + " (#elim " +
				 to_string(nnzs_thr-nnzs_aln) + " " +
				 to_string((static_cast<double>(nnzs_thr-nnzs_aln)/
							nnzs_thr)*100.0)
				 + "%)");
	// C.ParallelWriteMM(params.align_file, true, CkOutputHandler<NT>());	
}


template void pw_aln_batch<CommonKmerLoc>
(std::shared_ptr<DistFastaData>,
 combblas::SpParMat <uint64_t, CommonKmerLoc,
 combblas::SpDCCols <uint64_t, CommonKmerLoc>> &,
 PWAlign &, params_t &, int, int);

template void pw_aln_batch<CommonKmerLight>
(std::shared_ptr<DistFastaData>,
 combblas::SpParMat <uint64_t, CommonKmerLight,
 combblas::SpDCCols <uint64_t, CommonKmerLight>> &,
 PWAlign &, params_t &, int, int);

}
