// Created by Saliya Ekanayake on 2019-09-03.

#pragma once

#include <vector>

#include "seqan/score.h"
#include "seqan/sequence.h"

#include "PWAlign.hpp"



namespace
pastis
{



// uses location representation
class
SeqanXdropAligner : public PWAlign
{

private:

	std::vector<seqan::Peptide> rseqs_;
	std::vector<seqan::Peptide> cseqs_;
	seqan::Blosum62				blosum62_;
	ushort						seed_len_;
	int							xdrop_;
	int							seed_cnt_;



public:

	SeqanXdropAligner (int gap_open, int gap_ext,
					   ushort seed_len, int xdrop, int seed_cnt) :
		PWAlign(), blosum62_(gap_open, gap_ext),
		seed_len_(seed_len), xdrop_(xdrop), seed_cnt_(seed_cnt)
	{
	}



	~SeqanXdropAligner ()
	{
	}


	void
	construct_seqs (std::shared_ptr<DistFastaData> dfd) override;



	void
	construct_seqs_bl (std::shared_ptr<DistFastaData> dfd) override;



	void
	aln_batch (std::tuple<uint64_t, uint64_t, CommonKmerLoc *> *mattuples,
			   uint64_t beg, uint64_t end, 
			   uint64_t bl_roffset, uint64_t bl_coffset,
			   const params_t &params) override;
};


	
}
