// Created by Saliya Ekanayake on 2/7/19.

#pragma once

#include <sched.h>

#include <cstdint>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <utility>
#include <vector>

#include "CombBLAS/CombBLAS.h"

#include "Alphabet.hpp"
#include "ParallelOps.hpp"
#include "sr.hpp"

typedef unsigned short ushort;
typedef unsigned char uchar;
typedef std::vector<uint64_t> uvec_64;
typedef std::vector<uint64_t> uvec_16;


extern shared_ptr<pastis::ParallelOps> parops;



namespace
pastis
{



// pastis parameters
struct
params_t
{
	std::string	input_file;
	uint64_t	input_overlap;
	uint64_t	seq_count;
	int			gap_open;
	int			gap_ext;
	ushort		klength;
	ushort		kstride;

	bool		write_overlaps;
	std::string	overlap_file;
	bool		add_substitue_kmers;
	int			subk_count;

	// alignment
	enum PwAln
		{
			ALN_NONE = 0,
			ALN_SEQAN_FULL,
			ALN_SEQAN_XDROP,
			ALN_SEQAN_BANDED,
			ALN_ADEPT_GPUBSW
		};
	std::string	align_file;
	int			afreq;
	PwAln		pw_aln;
	int			aln_seqan_banded_hw;
	int			aln_seqan_xdrop;
	uint64_t	aln_batch_sz;
	double 		aln_cov_thr;	// alignment cov threshold
	int			aln_ani_thr;	// alignment ANI threshold

	std::string idx_map_file;

	Alphabet	*alph;
	std::string	 alph_str;

	int seed_count;

	int		ckthr;
	float	mosthr;

	int br, bc;					// blocked multiplication dimensions

	bool set_aff;
	std::vector<int> omp_thd_aff; // #cores per task
	std::vector<int> gpu_thd_aff; // #gpus per task

	bool pb;					// pre-blocking (i.e., overlap)
	enum LoadBal
		{
			LB_IDX = 0,
			LB_TRG
		};
	LoadBal lb;

	bool stats;
};



// parallel sparse matrix type shorthands
template <class NT>
struct
PSpMat
{
	typedef combblas::SpTuples <uint64_t, NT> Tuples;
	typedef combblas::SpDCCols <uint64_t, NT> DCCols;
	typedef combblas::SpParMat <uint64_t, NT, DCCols> MPI_DCCols;
	typedef std::tuple<uint64_t, uint64_t, NT *> ref_tuples;
};



// a nonzero element in the seq-by-kmer matrix
struct
MatrixEntry
{
	short cost;					// for score in computing kmer similarity
    ushort offset;				// the location of kmer in seq



	MatrixEntry () :
		cost(0), offset(0)
	{
	}



	MatrixEntry (short cost, ushort offset):
		cost(cost), offset(offset)
	{
	}



	friend
	std::ostream &
	operator<< (std::ostream& os, const MatrixEntry& me)
	{
		os << "(" << me.cost << ", " << me.offset << ")";
		return os;
	}



	MatrixEntry
	operator+ (const MatrixEntry &me) const
	{
		return me;
	}



	bool
	operator< (const MatrixEntry &me) const
	{
		return cost < me.cost;
	}	
};



// a nonzero in the output matrix - stores two matches with their locations in
// the seq pair
struct
CommonKmerLoc
{
	ushort count;
	short score;
	float score_aln;
	// store two overlaps per seq pair
	std::pair<ushort, ushort> first; // first kmer's locs in two seqs
    std::pair<ushort, ushort> second; // second kmer's locs in two seqs



	CommonKmerLoc () :
		count(1), score(-1), score_aln(0.0f)
	{
	}



	// CommonKmerLoc (short score, float score_aln) :
	// 	score(score),
	// 	score_aln(score_aln)
	// {
	// }



	explicit
	CommonKmerLoc (ushort count) :
		count(count)
	{
    }



	friend
	std::ostream &
	operator<< (std::ostream &os, const CommonKmerLoc &m)
	{
		os << "| " << m.count << " (" << m.first.first << "," << m.first.second
		   << ")(" <<
		   m.second.first << "," << m.second.second << ") | "
		   << m.score;
		return os;
    }



	operator float()
	{
		return score_aln;
	}
};



// a nonzero in the output matrix - does not keep locations of the kmers
struct
CommonKmerLight
{
	ushort count;				// overlap count
	short score;				// overlap score
	float score_aln;			// alignment score


	
	CommonKmerLight () :
		count(1), score(-1), score_aln(0.0f)
	{
	}



	// CommonKmerLight (short score, float score_aln) :
	// 	score(score),
	// 	score_aln(score_aln)
	// {
	// }



	explicit
	CommonKmerLight (ushort count) :
		count(count)
	{
    }



	friend
	std::ostream &
	operator<< (std::ostream &os, const CommonKmerLight &m)
	{
		os << m.count << " " << m.score << " " << m.score_aln;
		return os;
    }



	operator float()
	{
		return score_aln;
	}
};



template <typename NT>
struct
CkOutputHandler
{
    template <typename c, typename t>
	void save(std::basic_ostream<c,t> &os,
			  const NT &v, uint64_t row, uint64_t col)
	{
		os << " -> " << v.count << " " << v.score << " " << v.score_aln;
	}
};



// thread data for blocked multiplication
template <typename MatrixEntry,
		  typename T_OUT>
struct
MultData
{
	using SR = typename std::conditional<
			std::is_same<T_OUT, CommonKmerLight>::value,
			KmerIntersectLight<MatrixEntry, T_OUT>,
			KmerIntersectLoc<MatrixEntry, T_OUT>>::type;
	using DER_IN = typename combblas::SpDCCols<uint64_t, MatrixEntry>;
	using DER_OUT = typename combblas::SpDCCols<uint64_t, T_OUT>;


	combblas::BlockSpGEMM<uint64_t, MatrixEntry, DER_IN,
						  MatrixEntry, DER_IN> *bsp_;
	uint64_t roffset_, coffset_;
	int bri_, bci_;				// holds last computed block id
	int br_, bc_;
	combblas::SpParMat<uint64_t, T_OUT, DER_OUT> C_[2];
	int curC_;
	params_t *params_;

	std::vector<uint64_t> roffsets_;
	std::vector<uint64_t> coffsets_;
	bool partial_block_;
	


	MultData (combblas::BlockSpGEMM<uint64_t, MatrixEntry,
			  DER_IN, MatrixEntry, DER_IN> *bspgemm,
			  params_t *params
			  ) :
		bsp_(bspgemm), bri_(0), bci_(-1),
		br_(params->br), bc_(params->bc), curC_(1),
		params_(params)
	{
		roffsets_ = std::move(bspgemm->getBlockOffsets(true));
		coffsets_ = std::move(bspgemm->getBlockOffsets(false));
	}



	bool
	hasNext ()
	{
		return !((bri_ == br_-1) && (bci_ == bc_-1));
	}



	void
	computeNext (int nthds = -1)
	{
		getNextBlockIds();
	
		int numThreads = 1;
		#pragma omp parallel
    	{
			numThreads = omp_get_num_threads();
    	}
		
		if (nthds != -1)
			omp_set_num_threads(nthds);

		auto t_begin = std::chrono::system_clock::now();
		parops->logger->log("[computeNext] BEGIN bri " + to_string(bri_) +
							" bci " + to_string(bci_));

		parops->tp->start_timer("sparse");
		parops->tp->start_timer("mat-mult");

		curC_ = 1 - curC_;
		C_[curC_] = std::move(bsp_->template getBlockId<
			SR, T_OUT, combblas::SpDCCols<uint64_t, T_OUT > >
							  (bri_, bci_, roffset_, coffset_));

		parops->tp->stop_timer("sparse");
		parops->tp->stop_timer("mat-mult");

		parops->logger->log("[computeNext] END bri " + to_string(bri_) +
							" bci " + to_string(bci_));
		parops->logger->log("[computeNext] computed nnz " +
							to_string(C_[curC_].getnnz()));
		
		double t_block_pure =
			ms_t(std::chrono::system_clock::now()-t_begin).count();
		parops->info
			("Pure block took " + to_string(t_block_pure/1e3) + " seconds" +
			 " #threads used " +
			 to_string(nthds == -1 ? numThreads : nthds));

		if (nthds != -1)		// revert back to full
			omp_set_num_threads(numThreads);
	}



	// for triangular load balancing
	void
	computeNext_trg (int nthds = -1)
	{
		bool skip_block = false, full_block = false;
		do
		{
			getNextBlockIds();		
			skip_block = roffsets_[bri_] >= (coffsets_[bci_+1]-1);
			full_block = coffsets_[bci_] > (roffsets_[bri_+1]-1);
		} while (skip_block && hasNext());
		partial_block_ = !full_block && !skip_block;
		int numThreads = 1;
		#pragma omp parallel
    	{
			numThreads = omp_get_num_threads();
    	}
		
		if (nthds != -1)
			omp_set_num_threads(nthds);
		
		auto t_begin = std::chrono::system_clock::now();
		parops->logger->log("[computeNext_trg] BEGIN bri " + to_string(bri_) +
							" bci " + to_string(bci_));

		parops->tp->start_timer("sparse");
		parops->tp->start_timer("mat-mult");

		curC_ = 1 - curC_;
		C_[curC_] = std::move(bsp_->template getBlockId<
			SR, T_OUT, combblas::SpDCCols<uint64_t, T_OUT > >
							  (bri_, bci_, roffset_, coffset_));

		parops->tp->stop_timer("sparse");
		parops->tp->stop_timer("mat-mult");

		parops->logger->log("[computeNext_trg] END bri " + to_string(bri_) +
							" bci " + to_string(bci_));
		parops->logger->log("[computeNext_trg] computed nnz " +
							to_string(C_[curC_].getnnz()));
		
		double t_block_pure =
			ms_t(std::chrono::system_clock::now()-t_begin).count();
		parops->info
			("Pure block took " + to_string(t_block_pure/1e3) + " seconds" +
			 " #threads used " +
			 to_string(nthds == -1 ? numThreads : nthds));

		if (nthds != -1)		// revert back to full
			omp_set_num_threads(numThreads);
	}



	void
	getNextBlockIds (int &bri, int &bci)
	{
		bri = bri_;
		bci = bci_;
		if (bci_ == -1)			// first time
			bci = 0;
		else
		{
			if (bci_ == bc_-1)
			{
				bci = 0;
				bri = bri_ + 1;
			}
			else
				bci = bci_ + 1;
		}
		
	}
   



	// prevent calling computed block matrices' destructors?
	~MultData ()
	{
		
	}



private:

	void
	getNextBlockIds ()
	{
		if (bci_ == -1)			// first time
			bci_ = 0;
		else
		{
			if (bci_ == bc_-1)
			{
				bci_ = 0;
				++bri_;
			}
			else
				++bci_;
		}
	}
};
		  
	
}
