// Created by Saliya Ekanayake on 2019-07-05.


#include "../../include/pw/SeedExtendXdrop.hpp"


//template<typename TSequenceValue, typename TSpec>
//SeedExtendXdrop<TSequenceValue, TSpec>::SeedExtendXdrop(
SeedExtendXdrop::SeedExtendXdrop(
    seqan::Blosum62 scoring_scheme,
    seqan::Blosum62 scoring_scheme_simple,
    ushort seed_length, int xdrop, int seed_count):
    PairwiseFunction(),
    scoring_scheme(scoring_scheme),
    scoring_scheme_simple(scoring_scheme_simple),
    seed_length(seed_length), xdrop(xdrop), seed_count(seed_count){

}

//template<typename TSequenceValue, typename TSpec>
//void SeedExtendXdrop<TSequenceValue, TSpec>::apply(
void SeedExtendXdrop::apply(
    uint64_t l_col_idx, uint64_t g_col_idx,
    uint64_t l_row_idx, uint64_t g_row_idx,
    seqan::Peptide *seq_h, seqan::Peptide *seq_v,
    pastis::CommonKmers &cks, std::stringstream& ss) {

  AlignmentInfo ai[2];
  for (int count = 0; count < seed_count; ++count) {
    // row sequence is the same thing as vertical sequence
    ushort l_row_seed_start_offset = (count == 0) ? cks.first.first
                                                  : cks.second.first;
    // col sequence is the same thing as horizontal sequence
    ushort l_col_seed_start_offset = (count == 0) ? cks.first.second
                                                  : cks.second.second;

    // Seed creation params are:
    // horizontal seed start offset, vertical seed start offset, length
    TSeed seed(l_col_seed_start_offset, l_row_seed_start_offset,
               seed_length);
    auto start_time = std::chrono::system_clock::now();
    extendSeed(seed, *seq_h, *seq_v, seqan::EXTEND_BOTH, scoring_scheme_simple,
               xdrop,
               seqan::GappedXDrop());
    auto end_time = std::chrono::system_clock::now();
    add_time("XA:extend_seed", (ms_t(end_time - start_time)).count());

    seqan::Align<seqan::Peptide> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), infix(*seq_h, beginPositionH(seed),
                                      endPositionH(seed)));
    assignSource(row(align, 1), infix(*seq_v, beginPositionV(seed),
                                      endPositionV(seed)));

    /*! Note. This aligns the extended seeds globally, NOT the original
     * two sequences.
     *
     * It seems kind of a waste to have to do the alignment
     * again after xdrop seed extension but that's the only
     * way to get the alignment info in SeqAn.
     * See https://seqan.readthedocs.io/en/master/Tutorial/Algorithms/SeedExtension.html
     */
    start_time = std::chrono::system_clock::now();
    globalAlignment(align, scoring_scheme);
    end_time = std::chrono::system_clock::now();
    add_time("XA:global_alignment", (ms_t(end_time - start_time)).count());

    // Compute the statistics of the alignment.

    start_time = std::chrono::system_clock::now();
    computeAlignmentStats(ai[count].stats, align, scoring_scheme);
    end_time = std::chrono::system_clock::now();
    add_time("XA:compute_stats", (ms_t(end_time - start_time)).count());

    ai[count].seq_h_length = length(*seq_h);
    ai[count].seq_v_length = length(*seq_v);
    ai[count].seq_h_seed_length = static_cast<ushort>(seed._endPositionH -
                                                      seed._beginPositionH);
    ai[count].seq_v_seed_length = static_cast<ushort>(seed._endPositionV -
                                                      seed._beginPositionV);
    ai[count].seq_h_g_idx = g_col_idx;
    ai[count].seq_v_g_idx = g_row_idx;
  }

  /* Hard coding quality constraints for now */

  // TODO - Saliya
  // For now only keeps the largest alignment > 30% identity.
  // Incorporate length coverage restrictions later.

  AlignmentInfo max_ai = ai[0];
  if (seed_count > 2) {
    max_ai = ai[0].stats.alignmentIdentity > ai[1].stats.alignmentIdentity
    ? ai[0] : ai[1];
  }

//  if (max_ai.stats.alignmentIdentity >= 30.0){
//    alignments.push_back(max_ai);
//  }


  double alen_minus_gapopens = (max_ai.stats.alignmentLength - max_ai.stats.numGapOpens) * 1.0;
  ss << g_col_idx << "," << g_row_idx << "," << max_ai.stats.alignmentIdentity
     << "," << max_ai.seq_h_length << "," << max_ai.seq_v_length
     << "," << max_ai.seq_h_seed_length  << "," << max_ai.seq_v_seed_length
     << "," << max_ai.stats.numGapOpens
     << "," << alen_minus_gapopens / max_ai.seq_h_length
     << "," << alen_minus_gapopens / max_ai.seq_v_length
	 << "," << cks.count
	 << std::endl;
}



// @NOTE This is hard-coded to the number of seeds being <= 2
void
SeedExtendXdrop::apply_batch
(
    seqan::StringSet<seqan::Gaps<seqan::Peptide>> &seqsh,
	seqan::StringSet<seqan::Gaps<seqan::Peptide>> &seqsv,
	uint64_t								*lids,
	uint64_t								 col_offset,
	uint64_t								 row_offset,
	PSpMat<pastis::CommonKmers>::ref_tuples *mattuples,
	std::ofstream							&lfs,
	double									 thr_cov,
	int										 thr_ani
)
{
	seqan::ExecutionPolicy<seqan::Parallel, seqan::Vectorial> exec_policy;

	int numThreads = 1;
	#ifdef THREADED
	#pragma omp parallel
    {
      	numThreads = omp_get_num_threads();
    }
	#endif

	uint64_t npairs = seqan::length(seqsh);
	setNumThreads(exec_policy, numThreads);

	lfs << "processing batch of size " << npairs << " with "
		<< numThreads << " threads " << std::endl;

	// for multiple seeds we store the seed with the highest identity
	AlignmentInfo *ai = new AlignmentInfo[npairs];
	std::pair<ushort, ushort> *seedlens = new std::pair<ushort, ushort>[npairs];

	for (int count = 0; count < seed_count; ++count)
	{
		auto start_time = std::chrono::system_clock::now();

		seqan::StringSet<seqan::Gaps<seqan::Peptide>> seqsh_ex;
		seqan::StringSet<seqan::Gaps<seqan::Peptide>> seqsv_ex;
		resize(seqsh_ex, npairs, seqan::Exact{});
		resize(seqsv_ex, npairs, seqan::Exact{});
		
		// extend the current seed and form a new gaps object
		#pragma omp parallel for
		for (uint64_t i = 0; i < npairs; ++i)
		{
			pastis::CommonKmers *cks = std::get<2>(mattuples[lids[i]]);
			ushort l_row_seed_start_offset =
				(count == 0) ? cks->first.first : cks->second.first;
			ushort l_col_seed_start_offset =
				(count == 0) ? cks->first.second : cks->second.second;

			TSeed seed(l_col_seed_start_offset, l_row_seed_start_offset,
					   seed_length);
			extendSeed(seed, seqan::source(seqsh[i]), seqan::source(seqsv[i]),
					   seqan::EXTEND_BOTH, scoring_scheme_simple,
					   xdrop, seqan::GappedXDrop());
			assignSource(seqsh_ex[i],
						 infix(seqan::source(seqsh[i]),
							   beginPositionH(seed), endPositionH(seed)));
			assignSource(seqsv_ex[i],
						 infix(seqan::source(seqsv[i]),
							   beginPositionV(seed), endPositionV(seed)));
			seedlens[i].first = static_cast<ushort>(seed._endPositionH -
													seed._beginPositionH);
			seedlens[i].second = static_cast<ushort>(seed._endPositionV -
													 seed._beginPositionV);
		}

		auto end_time = std::chrono::system_clock::now();
    	add_time("XA:extend_seed", (ms_t(end_time - start_time)).count());

		
		start_time = std::chrono::system_clock::now();

		// alignment
		globalAlignment(exec_policy, seqsh_ex, seqsv_ex, scoring_scheme);
		
		end_time = std::chrono::system_clock::now();
    	add_time("XA:global_alignment", (ms_t(end_time - start_time)).count());


		start_time = std::chrono::system_clock::now();
		
		// stats
		if (count == 0)			// overwrite in the first seed
		{
			#pragma omp parallel for
			for (uint64_t i = 0; i < npairs; ++i)
			{
				computeAlignmentStats(ai[i].stats, seqsh_ex[i], seqsv_ex[i],
									  scoring_scheme);
				ai[i].seq_h_length = seqan::length(seqan::source(seqsh[i]));
				ai[i].seq_v_length = seqan::length(seqan::source(seqsv[i]));
				ai[i].seq_h_seed_length = seedlens[i].first;
				ai[i].seq_v_seed_length = seedlens[i].second;
				ai[i].seq_h_g_idx = col_offset + std::get<1>(mattuples[lids[i]]);
    			ai[i].seq_v_g_idx = row_offset + std::get<0>(mattuples[lids[i]]);
			}
		}
		else
		{
			#pragma omp parallel for
			for (uint64_t i = 0; i < npairs; ++i)
			{
				seqan::AlignmentStats stats;
				computeAlignmentStats(stats, seqsh_ex[i], seqsv_ex[i],
									  scoring_scheme);
				if (stats.alignmentIdentity > ai[i].stats.alignmentIdentity)
				{
					ai[i].stats				= stats;
					ai[i].seq_h_seed_length = seedlens[i].first;
					ai[i].seq_v_seed_length = seedlens[i].second;
				}
			}
		}

		end_time = std::chrono::system_clock::now();
    	add_time("XA:compute_stats", (ms_t(end_time - start_time)).count());
	}

	delete [] seedlens;

	auto start_time = std::chrono::system_clock::now();

	// stats
	#pragma omp parallel
	{
		#pragma omp for
		for (uint64_t i = 0; i < npairs; ++i)
		{
			seqan::AlignmentStats &stats = ai[i].stats;			
			double alen_minus_gapopens =
				stats.alignmentLength - stats.numGapOpens;

			// only keep alignments that meet coverage and ani criteria
			if (std::max((alen_minus_gapopens / ai[i].seq_h_length),
						 (alen_minus_gapopens / ai[i].seq_v_length)) >= thr_cov &&
				stats.alignmentIdentity >= thr_ani)
			{
				pastis::CommonKmers *cks = std::get<2>(mattuples[lids[i]]);
				cks->score_aln = (float)stats.alignmentIdentity / 100.0f;
				cks->score = 1;	// keep this
			}
		}
	}

	auto end_time = std::chrono::system_clock::now();
  	add_time("XA:string_op",
			 (ms_t(end_time - start_time)).count());

	delete [] ai;

	return;
}



// void
// SeedExtendXdrop::apply_batch_sc
// (
//     seqan::StringSet<seqan::Gaps<seqan::Peptide>> &seqsh,
// 	seqan::StringSet<seqan::Gaps<seqan::Peptide>> &seqsv,
// 	uint64_t *lids,
// 	uint64_t col_offset,
// 	uint64_t row_offset,
// 	PSpMat<pastis::CommonKmers>::Tuples &mattuples,
// 	std::ofstream &afs,
// 	std::ofstream &lfs
// )
// {
// 	seqan::ExecutionPolicy<seqan::Parallel, seqan::Vectorial> exec_policy;

// 	int numThreads = 1;
// 	#ifdef THREADED
// 	#pragma omp parallel
//     {
//       	numThreads = omp_get_num_threads();
//     }
// 	#endif

// 	uint64_t npairs = seqan::length(seqsh);
// 	setNumThreads(exec_policy, numThreads);

// 	lfs << "processing batch of size " << npairs << " with "
// 		<< numThreads << " threads " << std::endl;

// 	// for multiple seeds we store the seed with the highest identity
// 	seqan::String<int16_t> scores;
// 	seqan::String<int16_t> lens_h;
// 	seqan::String<int16_t> lens_v;
// 	resize(scores, npairs, seqan::Exact{});
// 	resize(lens_h, npairs, seqan::Exact{});
// 	resize(lens_v, npairs, seqan::Exact{});

// 	for (int count = 0; count < seed_count; ++count)
// 	{
// 		auto start_time = std::chrono::system_clock::now();

// 		seqan::StringSet<seqan::Peptide> seqsh_ex;
// 		seqan::StringSet<seqan::Peptide> seqsv_ex;
// 		resize(seqsh_ex, npairs, seqan::Exact{});
// 		resize(seqsv_ex, npairs, seqan::Exact{});
		
// 		// extend the current seed and form a new gaps object
// 		#pragma omp parallel for
// 		for (uint64_t i = 0; i < npairs; ++i)
// 		{
// 			pastis::CommonKmers &cks = mattuples.numvalue(lids[i]);
// 			ushort l_row_seed_start_offset =
// 				(count == 0) ? cks.first.first : cks.second.first;
// 			ushort l_col_seed_start_offset =
// 				(count == 0) ? cks.first.second : cks.second.second;

// 			TSeed seed(l_col_seed_start_offset, l_row_seed_start_offset,
// 					   seed_length);
// 			extendSeed(seed, seqsh[i], seqsv[i], seqan::EXTEND_BOTH,
// 					   scoring_scheme_simple, xdrop, seqan::GappedXDrop());
// 			seqsh_ex[i] = infix(seqsh[i], beginPositionH(seed),
// 								endPositionH(seed));
// 			seqsv_ex[i] = infix(seqsv[i], beginPositionV(seed),
// 								endPositionV(seed));
// 		}

// 		auto end_time = std::chrono::system_clock::now();
//     	add_time("XA:extend_seed", (ms_t(end_time - start_time)).count());

		
// 		start_time = std::chrono::system_clock::now();

// 		// alignment
// 		seqan::String<int16_t> scores_cur =
// 			globalAlignmentScore(exec_policy, seqsh_ex, seqsv_ex,
// 								 scoring_scheme);
		
// 		end_time = std::chrono::system_clock::now();
//     	add_time("XA:global_alignment", (ms_t(end_time - start_time)).count());


// 		start_time = std::chrono::system_clock::now();
		
// 		// stats
// 		if (count == 0)			// overwrite in the first seed
// 		{
// 			#pragma omp parallel for
// 			for (uint64_t i = 0; i < npairs; ++i)
// 			{
// 				scores[i] = scores_cur[i];
// 				lens_h[i] = seqan::length(seqsh_ex[i]);
// 				lens_v[i] = seqan::length(seqsv_ex[i]);
// 			}
// 		}
// 		else
// 		{
// 			#pragma omp parallel for
// 			for (uint64_t i = 0; i < npairs; ++i)
// 			{
// 				if (scores_cur[i] > scores[i])
// 				{
// 					scores[i] = scores_cur[i];
// 					lens_h[i] = seqan::length(seqsh_ex[i]);
// 					lens_v[i] = seqan::length(seqsv_ex[i]);
// 				}
// 			}
// 		}

// 		end_time = std::chrono::system_clock::now();
//     	add_time("XA:compute_stats", (ms_t(end_time - start_time)).count());
// 	}
	

// 	auto start_time = std::chrono::system_clock::now();

// 	// stats dump
// 	#pragma omp parallel
// 	{
// 		std::stringstream ss;

// 		#pragma omp for
// 		for (uint64_t i = 0; i < npairs; ++i)
// 		{
// 			if (scores[i] > 0)
// 			{
// 				int len_seqh = seqan::length(seqsh[i]);
// 				int len_seqv = seqan::length(seqsv[i]);
// 				ss << (col_offset + mattuples.colindex(lids[i])) << ","
// 				   << (row_offset + mattuples.rowindex(lids[i]))  << ","
// 				   << scores[i] << ","
// 				   << len_seqh << ","
// 				   << len_seqv << ","
// 				   << ((double)scores[i]/(double)len_seqh) << ","
// 				   << ((double)scores[i]/(double)len_seqv) << ","
// 				   << ((double)scores[i]/(double)lens_h[i]) << ","
// 				   << ((double)scores[i]/(double)lens_v[i])
// 				   << "\n";
// 			}
// 		}

// 		#pragma omp critical
// 		{
// 			afs << ss.str();
// 			afs.flush();
// 		}
// 	}

// 	auto end_time = std::chrono::system_clock::now();
//   	add_time("XA:string_op",
// 			 (ms_t(end_time - start_time)).count());

// 	return;
// }
