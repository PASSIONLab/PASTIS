// Created by Saliya Ekanayake on 10/15/19.

#include <unistd.h>

#include <algorithm>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#include "../../inc/kmer/KmerOps.hpp"
#include "../../inc/macros.hpp"
#include "../../inc/ScoreMat.hpp"
#include "../../inc/util.hpp"

using std::unordered_set;	using std::vector;	using std::string;
using std::accumulate;

extern shared_ptr<pastis::ParallelOps> parops;



namespace
pastis
{

static const char L2U_MASK = 0xDF;



static
uint64_t
add_kmers
(
    const char						*seq,
	ushort							 len,
	uint64_t						 start_offset,
	uint64_t						 end_offset_inclusive,
	ushort							 k,
	ushort							 s,
	Alphabet						&alp,
	uvec_64							&lcol_ids,
	vector<MatrixEntry>				&lvals,
	const shared_ptr<ParallelOps>	&parops,
	bool 							*kmer_universe,
	uint64_t						&nel_kmers
)
{	
	Blosum62 bsm62;
	auto num_kmers = static_cast<ushort>((floor((len-k)*1.0/s))+1);
	ushort base = alp.size;
	uint64_t kcode = 0;
	uint64_t count = 0;
	char cap_c;	
	uint64_t kmer_universe_size = pow(alp.size, k);
	char kmer_str[k+1];

	kmer_str[k] = '\0';
	for (uint64_t i = start_offset; i <= end_offset_inclusive-k+1; i += s)
	{
		kcode = 0;
		for (uint64_t j = i; j < i+k; ++j)
		{
			cap_c = *(seq + j) & L2U_MASK;	// lower to upper case
			kmer_str[j-i] = alp.al_map[cap_c];
			kcode		  = kcode*base + alp.char_to_code[cap_c];

			// encountered a char not in protein alphabet. yes it happens
			if (alp.char_to_code[cap_c] == Alphabet::capacity)
			{
				++nel_kmers;
				goto skip_kmer;
			}
		}
		
		lcol_ids.push_back(kcode);
		kmer_universe[kcode] = true;
		lvals.emplace_back(bsm62.self_score(kmer_str),
						   static_cast<ushort>(i-start_offset));
		++count;
		
	skip_kmer:
		;
	}

	return count;
}

	


PSpMat<MatrixEntry>::MPI_DCCols
generate_A
(
    uint64_t						 seq_count,
	shared_ptr<DistFastaData>		&dfd,
	ushort							 k,
	ushort							 s,
	Alphabet						&alph,
	unordered_set<Kmer, Kmer>		&local_kmers
)
{	
	char				*buff;
    ushort				 len;
    uint64_t			 start_offset, end_offset_inclusive;
    uvec_64				 lrow_ids, lcol_ids;
	vector<MatrixEntry>	 lvals;
	uint64_t			 offset = dfd->global_start_idx();
    FastaData			*lfd	= dfd->lfd();
	uint64_t			 nel_kmers = 0;

	parops->logger->log("Generating seq-by-kmer matrix.");
	parops->logger->log("Parsing kmers.");

	// @OGUZ-TODO need reduction on when greater than MAX_INT
	uint64_t tmp = pow(alph.size, k);
	assert (tmp >= 0 && tmp <= std::numeric_limits<int>::max() &&
			"kmer universe size greater than MAX_INT. "
			"Unable to do the MPI reduction. Abort.");
	int kmer_universe_size = tmp;
	bool *kmer_universe = new bool[kmer_universe_size]();

	#if PASTIS_DBG_LVL > 0
	parops->bytes_alloc += kmer_universe_size * sizeof(bool);
	#endif

	parops->logger->log("approximate memory in usage " +
						gb_str(parops->bytes_alloc));

	parops->tp->start_timer("kmerop:gen_A:loop_add_kmers");
	for (uint64_t lseq_idx = 0; lseq_idx < lfd->local_count(); ++lseq_idx)
	{
		// parops->logger->debug("seq id " + to_string(lseq_idx));
		buff = lfd->get_sequence(lseq_idx, len, start_offset,
								 end_offset_inclusive);
		auto num_kmers =
			add_kmers(buff, len, start_offset, end_offset_inclusive, k, s,
					  alph, lcol_ids, lvals, parops, kmer_universe, nel_kmers);
		lrow_ids.insert(lrow_ids.end(), num_kmers, lseq_idx + offset);
	}
	parops->tp->stop_timer("kmerop:gen_A:loop_add_kmers");

	parops->logger->log("#eliminated kmers due to foreign chars " +
						to_string(nel_kmers));
	parops->logger->log("kmer universe size " + to_string(kmer_universe_size));
	parops->logger->log("kmer_universe arr size " +
						mb_str(sizeof(bool)*kmer_universe_size));
	parops->logger->log("#elems in lrow_ids lcol_ids lvals: " +
						to_string(lrow_ids.size()) +
						" size each (lrow_ids lcol_ids) " +
						mb_str(sizeof(uint64_t)*lrow_ids.size()) + 
						" size lvals " +
						mb_str(sizeof(MatrixEntry)*lvals.size()));
	#if PASTIS_DBG_LVL > 0
	parops->bytes_alloc += lrow_ids.capacity() * sizeof(uint64_t);
	parops->bytes_alloc += lcol_ids.capacity() * sizeof(uint64_t);
	parops->bytes_alloc += lvals.capacity() * sizeof(MatrixEntry);
	#endif

	parops->logger->log("approximate memory in usage " +
						gb_str(parops->bytes_alloc));

	// @OGUZ-TODO Allreduce with size larger than INT_MAX
	MPI_Allreduce(MPI_IN_PLACE,
				  kmer_universe,
				  kmer_universe_size,
				  MPI_CXX_BOOL,
				  MPI_LOR,
				  MPI_COMM_WORLD);	

	uint64_t unique_kmer_count = 0;
	unique_kmer_count = accumulate(kmer_universe,
								   kmer_universe+kmer_universe_size, 0);

	parops->logger->log("unique kmer count " + to_string(unique_kmer_count));

	uint64_t q = unique_kmer_count / parops->g_np;
	uint64_t r = unique_kmer_count - (q*parops->g_np);
	uint64_t skip_count =
		(parops->g_rank < r) ? ((q+1)*parops->g_rank) : (q*parops->g_rank+r);
	uint64_t my_count = (parops->g_rank < r) ? q+1 : q;
	uint64_t my_end = my_count + skip_count;
	uint64_t count = 0;

	parops->logger->log("computing my portion of kmers");
	
	// @OGUZ-REVIEW parallelize this kmer_universe_size might get quite large
	for (size_t i = 0; i < kmer_universe_size; ++i)
	{
		if (kmer_universe[i])
			++count;
		else
			continue;

		if (count <= skip_count) // not mine
			continue;

		if (count <= my_end)
			local_kmers.emplace(i, k, alph);
		else					// cant be mine any more
			break;
	}

	delete[] kmer_universe;

	#if PASTIS_DBG_LVL > 0
	parops->bytes_alloc -= kmer_universe_size * sizeof(bool);
	#endif

	// form the seq-by-kmer matrix
	parops->logger->log("forming the seq-by-kmer matrix");
	assert(lrow_ids.size() == lcol_ids.size()
		   && lcol_ids.size() == lvals.size());

	uint64_t n_rows = seq_count;
	uint64_t n_cols = static_cast<uint64_t>(pow(alph.size, k));
	combblas::FullyDistVec<uint64_t, uint64_t> drows(lrow_ids, parops->grid);
	combblas::FullyDistVec<uint64_t, uint64_t> dcols(lcol_ids, parops->grid);
	combblas::FullyDistVec<uint64_t, MatrixEntry> dvals(lvals, parops->grid);
	
	parops->tp->start_timer("kmerop:gen_A:spMatA");
	PSpMat<MatrixEntry>::MPI_DCCols A(n_rows, n_cols,
									  drows, dcols, dvals, false);
	parops->tp->stop_timer("kmerop:gen_A:spMatA");

	#if PASTIS_DBG_LVL > 0
	parops->bytes_alloc += A.getlocalnnz() * sizeof(uint64_t);
	parops->bytes_alloc += A.getlocalnnz() * sizeof(MatrixEntry);
	parops->bytes_alloc += A.getlocalcols() * sizeof(uint64_t) * 2;
	#endif

	parops->logger->log("approximate memory in usage " +
						gb_str(parops->bytes_alloc));


	return A;
	// return PSpMat<MatrixEntry>::MPI_DCCols();
}

}
