// Created by Saliya Ekanayake on 2/7/19.

#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "CombBLAS/CombBLAS.h"

#include "Alphabet.hpp"

typedef unsigned short ushort;
typedef unsigned char uchar;
typedef std::vector<uint64_t> uvec_64;
typedef std::vector<uint64_t> uvec_16;



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

	
}
