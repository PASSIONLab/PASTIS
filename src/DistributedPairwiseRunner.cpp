//
// Created by Saliya Ekanayake on 2019-02-19.
// Modified by Aydin Buluc on 2019-12-29 
//


#include "../include/DistributedPairwiseRunner.hpp"
#include <atomic>         // std::atomic, std::atomic_flag, ATOMIC_FLAG_INIT


DistributedPairwiseRunner::DistributedPairwiseRunner(
    const std::shared_ptr<DistributedFastaData> dfd,
    PSpMat<pastis::CommonKmers>::DCCols * localmat,
	PSpMat<pastis::CommonKmers>::MPI_DCCols *glmat,
    int afreq,
    uint64_t rowoffset, uint64_t coloffset,
    const std::shared_ptr<ParallelOps> &parops)
    : dfd(dfd), gmat(glmat), spSeq(localmat), row_offset(rowoffset),
	  col_offset(coloffset), afreq(afreq), parops(parops) {

}


void DistributedPairwiseRunner::write_overlaps(const char *file) {

  uint64_t local_nnz_count = 0;
  uint64_t local_top_triangle_count = 0;
  std::stringstream ss;

  std::ofstream afs;
  afs.open(file);

  if (parops->world_proc_rank == 0) {
    ss << "g_col_idx,g_row_idx,common_kmer_count" << std::endl;
  }
  ushort l_max_common_kmers = 0;
  for (auto colit = spSeq->begcol(); colit != spSeq->endcol(); ++colit) {
    // iterate over columns
    auto l_col_idx = colit.colid(); // local numbering
    uint64_t g_col_idx = l_col_idx + col_offset;

    for (auto nzit = spSeq->begnz(colit); nzit < spSeq->endnz(colit); ++nzit) {
      auto l_row_idx = nzit.rowid();
      uint64_t g_row_idx = l_row_idx + row_offset;

      ++local_nnz_count;

      /*!
       * Note. the cells means the process grid cells.
       * We only want to compute the top triangles of any grid cell.
       * Further, we want cell diagonals only for cells that are on the
       * top half of the grid excluding the grid's main diagonal cells
       */
      if (l_col_idx < l_row_idx){
        continue;
      }

      if (l_col_idx == l_row_idx && g_col_idx <= g_row_idx){
        continue;
      }

      pastis::CommonKmers cks = nzit.value();
      if (cks.count > l_max_common_kmers){
        l_max_common_kmers = cks.count;
      }

      ++local_top_triangle_count;
      // ss << g_col_idx << "," << g_row_idx << "," << cks.count << "\n";
	  afs << g_row_idx << " " << g_col_idx << "\n";
    }
  }

  afs.close();

  ushort g_max_common_kmers = 0;
  MPI_Reduce(&l_max_common_kmers, &g_max_common_kmers, 1,
             MPI_UINT16_T, MPI_MAX, 0, MPI_COMM_WORLD);
  if (parops->world_proc_rank == 0){
    std::printf("  Max common kmers %d\n", g_max_common_kmers);
  }
  // std::string overlaps_str = ss.str();
  // parops->write_file_in_parallel(file, overlaps_str);
}

void DistributedPairwiseRunner::run(PairwiseFunction *pf, const char* file, std::ofstream& lfs, int log_freq) {
  /*! There are two types of rows and columns below.
   * The sequences are arranged as an NxN matrix in
   * mat (this is not how it's stored internally).
   * This NxN is distributed over a grid of
   * sqrt(P) x sqrt (P), where P is the total number
   * of processes. Anything to do with the grid will
   * be prefixed by gr_*/


  std::ofstream af_stream;
  af_stream.open(file);

  uint64_t local_nnz_count = spSeq->getnnz();
  std::atomic<uint64_t> current_nnz_count(0);

  lfs << "Local nnz count: " << local_nnz_count << std::endl;

  int numThreads = 1;	// default case
#ifdef THREADED
#pragma omp parallel
    {
      numThreads = omp_get_num_threads();
    }
#endif

  std::vector<std::stringstream> ss(numThreads);
  if(parops->world_proc_rank == 0)
    af_stream << "g_col_idx,g_row_idx,pid,col_seq_len,row_seq_len,"
		"col_seq_align_len,row_seq_align_len, num_gap_opens, "
		"col_seq_len_coverage, row_seq_len_coverage,common_count" << std::endl;

  std::atomic<uint64_t> line_count(0);
  uint64_t nalignments = 0;
  PSpMat<pastis::CommonKmers>::Tuples mattuples(*spSeq);
  
  #pragma omp parallel for reduction(+:nalignments)
  for(uint64_t i=0; i< local_nnz_count; i++)
  {
	auto l_row_idx = mattuples.rowindex(i);
	auto l_col_idx = mattuples.colindex(i);
	uint64_t g_col_idx = l_col_idx + col_offset;
	uint64_t g_row_idx = l_row_idx + row_offset;		  

	seqan::Peptide *seq_h = dfd->col_seq(l_col_idx);  
	seqan::Peptide *seq_v = dfd->row_seq(l_row_idx);

	current_nnz_count++;
	if (current_nnz_count % log_freq == 0){
		#pragma omp critical
		{
			  auto t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
			  lfs << "  (" << current_nnz_count << "/" << local_nnz_count << ") -- "
			  << std::setprecision(2) << (1.0*current_nnz_count / local_nnz_count)
			  << "% done. " << std::ctime(&t);
			  lfs.flush();
		}
	}	

	/*!
	 * Note. the cells means the process grid cells.
	 * We only want to compute the top triangles of any grid cell.
	 * Further, we want cell diagonals only for cells that are on the
	 * top half of the grid excluding the grid's main diagonal cells
	 * */
	if (l_col_idx < l_row_idx){
		continue;
	}
	if (l_col_idx == l_row_idx && g_col_idx <= g_row_idx){
		continue;
	}

	pastis::CommonKmers cks = mattuples.numvalue(i);


	int myThread = 0;
#ifdef THREADED
	myThread = omp_get_thread_num();
#endif

	++nalignments;
	pf->apply(l_col_idx, g_col_idx, l_row_idx, g_row_idx, seq_h, seq_v, cks, ss[myThread]);
	line_count++;

	if (line_count % afreq == 0)
	{
		#pragma omp critical
		{
		  af_stream << ss[myThread].str();
		  af_stream.flush();
		  ss[myThread].str(std::string());
		}
	}
  }

  pf->nalignments = nalignments;

  lfs << "  (" << current_nnz_count << "/" << local_nnz_count << ") -- "
      << "100% done." << std::endl;
  lfs << "#alignments run " << nalignments << std::endl;

  pf->print_avg_times(parops, lfs);
  
//  if(parops->world_proc_rank == 0) {
//  }
//  ushort g_max_common_kmers = 0;
//  MPI_Reduce(&l_max_common_kmers, &g_max_common_kmers, 1,
//             MPI_UINT16_T, MPI_MAX, 0, MPI_COMM_WORLD);
//  if (parops->world_proc_rank == 0){
//    std::printf("  Max common kmers %d\n", g_max_common_kmers);
//  }
//  std::string align_str = ss.str();
//  parops->write_file_in_parallel(file, align_str);

  for(int i=0; i< numThreads; ++i)
  {
  	af_stream << ss[i].str();
  }
  af_stream.flush();
  af_stream.close();
}



void
DistributedPairwiseRunner::run_batch
(
    PairwiseFunction	*pf,
	const std::string	&aln_file,
	std::ofstream&		 lfs,
	int					 log_freq,
	int					 ckthr,
	float				 mosthr,
	TraceUtils 			 tu,
	bool				 score_only
)
{
	uint64_t	local_nnz_count = spSeq->getnnz();
	int			batch_size		= 1e9;
	int			batch_cnt		= (local_nnz_count / batch_size) + 1;
	int			batch_idx		= 0;
	uint64_t	nalignments		= 0;

	// @TODO threaded
	PSpMat<pastis::CommonKmers>::ref_tuples *mattuples =
		new PSpMat<pastis::CommonKmers>::ref_tuples[local_nnz_count];
	uint64_t k = 0;
	auto dcsc = spSeq->GetDCSC();
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

	assert (k == local_nnz_count);
		
	lfs << "Local nnz count: " << local_nnz_count << std::endl;

	int numThreads = 1;
	#ifdef THREADED
	#pragma omp parallel
    {
      	numThreads = omp_get_num_threads();
    }
	#endif

	uint64_t *algn_cnts = new uint64_t[numThreads + 1];
	uint64_t nelims_ckthr = 0, nelims_mosthr = 0, nelims_both = 0;
	while (batch_idx < batch_cnt)
	{
		uint64_t beg = batch_idx * batch_size;
		uint64_t end = ((batch_idx + 1) * batch_size > local_nnz_count) ?
			local_nnz_count : ((batch_idx + 1) * batch_size);

		tu.print_str("batch idx " + std::to_string(batch_idx) + "/" +
					 std::to_string(batch_cnt) + " [" +
					 std::to_string(beg) + ", " +
					 std::to_string(end) + ")\n");

		memset(algn_cnts, 0, sizeof(*algn_cnts) * (numThreads + 1));

		uint64_t nelims_ckthr_cur = 0, nelims_mosthr_cur = 0,
			nelims_both_cur = 0;
		
		// count number of alignments in this batch
		#pragma omp parallel reduction(+:nelims_ckthr_cur,\
									   nelims_mosthr_cur,nelims_both_cur)
		{
			int tid = 0;
			#ifdef THREADED
			tid = omp_get_thread_num();
			#endif

			uint64_t algn_cnt = 0;

			#pragma omp for schedule(static, 1000)
			for (uint64_t i = beg; i < end; ++i)
			{
				// auto				l_row_idx = mattuples.rowindex(i);
				// auto				l_col_idx = mattuples.colindex(i);
				// uint64_t			g_col_idx = l_col_idx + col_offset;
				// uint64_t			g_row_idx = l_row_idx + row_offset;
				// pastis::CommonKmers cks		  = mattuples.numvalue(i);

				auto				 l_row_idx = std::get<0>(mattuples[i]);
				auto				 l_col_idx = std::get<1>(mattuples[i]);
				uint64_t			 g_col_idx = l_col_idx + col_offset;
				uint64_t			 g_row_idx = l_row_idx + row_offset;
				pastis::CommonKmers *cks	   = std::get<2>(mattuples[i]);
				
				if ((cks->count > ckthr) &&
					(cks->score > mosthr) &&
					(l_col_idx >= l_row_idx) &&
					(l_col_idx != l_row_idx || g_col_idx > g_row_idx))
					++algn_cnt;

				// stats purposes
				if ((l_col_idx >= l_row_idx) &&
					(l_col_idx != l_row_idx || g_col_idx > g_row_idx))
				{
					if (cks->count <= ckthr)
						++nelims_ckthr_cur;
					if (cks->score <= mosthr)
						++nelims_mosthr_cur;
					if (cks->count <= ckthr && cks->score <= mosthr)
						++nelims_both_cur;
				}
			}

			algn_cnts[tid + 1] = algn_cnt;
		}

		nelims_ckthr  += nelims_ckthr_cur;
		nelims_mosthr += nelims_mosthr_cur;
		nelims_both	  += nelims_both_cur;		

		// for (int i = 1; i < numThreads + 1; ++i)
		// 	lfs << "thread " << (i - 1) << ": " << algn_cnts[i] << " - ";
		// lfs << "\n" << "cur tot " << algn_cnts[numThreads] << " cum tot "
		// 	<< nalignments << std::endl;
		
		for (int i = 1; i < numThreads + 1; ++i)
			algn_cnts[i] += algn_cnts[i - 1];
		nalignments += algn_cnts[numThreads];

		if (algn_cnts[numThreads] == 0)
		{
			++batch_idx;
			continue;
		}
		
		// allocate StringSet
		seqan::StringSet<seqan::Gaps<seqan::Peptide>> seqsh;
		seqan::StringSet<seqan::Gaps<seqan::Peptide>> seqsv;
		resize(seqsh, algn_cnts[numThreads], seqan::Exact{});
		resize(seqsv, algn_cnts[numThreads], seqan::Exact{});
		uint64_t *lids = new uint64_t[algn_cnts[numThreads]];
		
		// fill StringSet
		#pragma omp parallel
		{
			int tid = 0;
			#ifdef THREADED
			tid = omp_get_thread_num();
			#endif

			uint64_t algn_idx = algn_cnts[tid];

			#pragma omp for schedule(static, 1000)
			for (uint64_t i = beg; i < end; ++i)
			{
				// auto				l_row_idx = mattuples.rowindex(i);
				// auto				l_col_idx = mattuples.colindex(i);
				// uint64_t			g_col_idx = l_col_idx + col_offset;
				// uint64_t			g_row_idx = l_row_idx + row_offset;
				// pastis::CommonKmers cks		  = mattuples.numvalue(i);

				auto				 l_row_idx = std::get<0>(mattuples[i]);
				auto				 l_col_idx = std::get<1>(mattuples[i]);
				uint64_t			 g_col_idx = l_col_idx + col_offset;
				uint64_t			 g_row_idx = l_row_idx + row_offset;
				pastis::CommonKmers *cks	   = std::get<2>(mattuples[i]);
				
				if ((cks->count > ckthr) &&
					(cks->score > mosthr) &&
					(l_col_idx >= l_row_idx) &&
					(l_col_idx != l_row_idx || g_col_idx > g_row_idx))
				{
					seqsh[algn_idx] =
						seqan::Gaps<seqan::Peptide>(*(dfd->col_seq(l_col_idx)));
					seqsv[algn_idx] =
						seqan::Gaps<seqan::Peptide>(*(dfd->row_seq(l_row_idx)));
					lids[algn_idx] = i;
					++algn_idx;
				}

				cks->score = 0; // will use it for pruning purposes
			}
		}
		

		// call aligner
		lfs << "calling aligner for batch idx " << batch_idx
			<< " cur #algnments " << algn_cnts[numThreads]
			<< " overall " << nalignments
			<< std::endl;
		// if (score_only)
		// 	pf->apply_batch_sc(seqsh, seqsv, lids, col_offset, row_offset,
		// 					   mattuples, af_stream, lfs);
		// else
		pf->apply_batch(seqsh, seqsv, lids, col_offset, row_offset,
						mattuples, lfs);
		
		
		delete [] lids;
		++batch_idx;
	}

	pf->nalignments = nalignments;
	pf->print_avg_times(parops, lfs);

	lfs << "#alignments run " << nalignments << std::endl;

	// Stats
	uint64_t nelims_ckthr_tot = 0, nelims_mosthr_tot = 0, nelims_both_tot = 0,
		nalignments_tot = 0;
	MPI_Reduce(&nelims_ckthr, &nelims_ckthr_tot, 1, MPI_UINT64_T,
			   MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&nelims_mosthr, &nelims_mosthr_tot, 1, MPI_UINT64_T,
			   MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&nelims_both, &nelims_both_tot, 1, MPI_UINT64_T,
			   MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&nalignments, &nalignments_tot, 1, MPI_UINT64_T,
			   MPI_SUM, 0, MPI_COMM_WORLD);
	tu.print_str("total nnzs in the output matrix " +
				 std::to_string(gmat->getnnz()) +
				 "\ntotal nnzs in strictly lower (or upper) mat " +
				 std::to_string((gmat->getnnz()-gmat->getncol())/2) + 
				 "\n  total alignments run " + std::to_string(nalignments_tot) +
				 "\n  eliminated due to common k-mer threshold " +
				 std::to_string(nelims_ckthr_tot) +
				 "\n  eliminated due to max overlap score threshold " +
				 std::to_string(nelims_mosthr_tot) +
				 "\n  eliminated due to both (intersection) " +
				 std::to_string(nelims_both_tot) +
				 "\n  eliminated overall with both (union) " +
				 std::to_string(nelims_ckthr_tot+nelims_mosthr_tot-
								nelims_both_tot) + "\n");

	// prune pairs that do not meet coverage criteria
	// std::string outfile = aln_file + std::string("-pruned-C.mtx");
	auto elim_cov = [] (pastis::CommonKmers &ck)
		{return ck.score == 0;};
	gmat->Prune(elim_cov);
	gmat->ParallelWriteMM(aln_file, true, pastis::CkOutputHandler());
	tu.print_str("nnzs in the pruned matrix " +
				 std::to_string(gmat->getnnz()) + "\n");
	
	delete [] algn_cnts;
	delete [] mattuples;

	return;
}
