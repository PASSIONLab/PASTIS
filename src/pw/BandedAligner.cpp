// Created by Saliya Ekanayake on 2019-09-03.

#include "../../include/pw/BandedAligner.hpp"

BandedAligner::BandedAligner(seqan::Blosum62 scoring_scheme,
    int banded_half_width) :

    PairwiseFunction(),
    scoring_scheme(scoring_scheme),
    banded_half_width(banded_half_width){

}

//template<typename TSequenceValue, typename TSpec>
//void SeedExtendXdrop<TSequenceValue, TSpec>::apply(
void BandedAligner::apply(
    uint64_t l_col_idx, uint64_t g_col_idx,
    uint64_t l_row_idx, uint64_t g_row_idx,
    seqan::Peptide *seq_h, seqan::Peptide *seq_v,
    CommonKmers &cks, std::stringstream& ss) {

  seqan::Align<seqan::Peptide> align;
  resize(rows(align), 2);
  assignSource(row(align, 0), *seq_h);
  assignSource(row(align, 1), *seq_v);

  auto start_time = std::chrono::system_clock::now();
  int score = localAlignment(align, scoring_scheme, -banded_half_width, banded_half_width);
  auto end_time = std::chrono::system_clock::now();
  add_time("BA:local_alignment", (ms_t(end_time - start_time)).count());

  AlignmentInfo ai;
  start_time = std::chrono::system_clock::now();
  computeAlignmentStats(ai.stats, align, scoring_scheme);
  end_time = std::chrono::system_clock::now();
  add_time("BA:compute_stats", (ms_t(end_time - start_time)).count());

  ai.seq_h_length = length(*seq_h);
  ai.seq_v_length = length(*seq_v);
  ai.seq_h_seed_length = (clippedEndPosition(row(align, 0)) - 1) - clippedBeginPosition(row(align, 0));
  ai.seq_v_seed_length = (clippedEndPosition(row(align, 1)) - 1) - clippedBeginPosition(row(align, 1));
  ai.seq_h_g_idx = g_col_idx;
  ai.seq_v_g_idx = g_row_idx;

  /* Hard coding quality constraints for now */

  // TODO - Saliya
  // For now only keeps the largest alignment > 30% identity.
  // Incorporate length coverage restrictions later.

//  if (max_ai.stats.alignmentIdentity >= 30.0){
  alignments.push_back(ai);
//  }

  ss << g_col_idx << "," << g_row_idx << "," << ai.stats.alignmentIdentity
     << "," << ai.seq_h_length << "," << ai.seq_v_length
     << "," << ai.seq_h_seed_length  << "," << ai.seq_v_seed_length << std::endl;
}