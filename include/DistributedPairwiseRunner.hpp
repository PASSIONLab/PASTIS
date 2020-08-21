//
// Created by Saliya Ekanayake on 2019-02-19.
// Modified by Aydin Buluc on 2019-12-29
//

#ifndef LBL_DAL_ALIGNER_H
#define LBL_DAL_ALIGNER_H

#include <seqan/align.h>
#include "Utils.hpp"
#include "DistributedFastaData.hpp"
#include "pw/PairwiseFunction.hpp"
#include "AlignmentInfo.hpp"
#include "kmer/CommonKmers.hpp"

//struct AlignmentInfo{
//  seqan::AlignmentStats stats;
//  ushort seq_h_length;
//  ushort seq_v_length;
//  ushort seq_h_seed_length;
//  ushort seq_v_seed_length;
//  uint64_t seq_h_g_idx;
//  uint64_t seq_v_g_idx;
//};

class DistributedPairwiseRunner {
public:
  DistributedPairwiseRunner(std::shared_ptr<DistributedFastaData> dfd,
                     PSpMat<pastis::CommonKmers>::DCCols * localmat,
							PSpMat<pastis::CommonKmers>::MPI_DCCols *glmat,
							int afreq,
		     uint64_t rowoffset, uint64_t coloffset,
                     const std::shared_ptr<ParallelOps> &parops);

//  uint64_t align_seqs();
  void write_overlaps(const char *file);
  void run(PairwiseFunction *pf, const char* file, std::ofstream& lfs, int log_freq);
  void run_batch(PairwiseFunction *pf, const char* file, std::ofstream& lfs,
				 int log_freq, int ckthr, float mosthr, TraceUtils tu,
				 bool score_only = false);

private:
  PSpMat<pastis::CommonKmers>::DCCols * spSeq;
  PSpMat<pastis::CommonKmers>::MPI_DCCols * gmat;
  uint64_t row_offset;  // local to global row id offset  
  uint64_t col_offset;	// ditto
  std::shared_ptr<DistributedFastaData> dfd;
  int afreq;
  std::shared_ptr<ParallelOps> parops;
};


#endif //LBL_DAL_ALIGNER_H
