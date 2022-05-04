// Created by Saliya Ekanayake on 2019-07-05.

#pragma once

#include <memory>
#include <tuple>

#include "../DistFastaData.hpp"
#include "../Types.hpp"



namespace
pastis
{




class
PWAlign
{

public:

	virtual
	~PWAlign ()
	{
	}

	

	virtual
	void
	construct_seqs (std::shared_ptr<DistFastaData> dfd) = 0;



	virtual
	void
	construct_seqs_bl (std::shared_ptr<DistFastaData> dfd) = 0;



	virtual
	void
	aln_batch (std::tuple<uint64_t, uint64_t, CommonKmerLight *> *mattuples,
			   uint64_t beg, uint64_t end, 
			   uint64_t bl_roffset, uint64_t bl_coffset,
			   const params_t &params)
	{
	}



	virtual
	void
	aln_batch (std::tuple<uint64_t, uint64_t, CommonKmerLoc *> *mattuples,
			   uint64_t beg, uint64_t end, 
			   uint64_t bl_roffset, uint64_t bl_coffset,
			   const params_t &params)
	{
	}



	virtual
	void
	aln_batch_ovlp (std::tuple<uint64_t, uint64_t,
					CommonKmerLight *> *mattuples,
					uint64_t beg, uint64_t end,
					uint64_t bl_roffset, uint64_t bl_coffset,
					MultData<MatrixEntry, CommonKmerLight> &md,
					const params_t &params)
	{
	}



	virtual
	void
	aln_batch_ovlp (std::tuple<uint64_t, uint64_t,
					CommonKmerLoc *> *mattuples,
					uint64_t beg, uint64_t end,
					uint64_t bl_roffset, uint64_t bl_coffset,
					MultData<MatrixEntry, CommonKmerLoc> &md,
					const params_t &params)
	{
	}




	virtual
	size_t
	rseq_len (uint64_t i)
	{
		return 0;
	}



	virtual
	size_t
	cseq_len (uint64_t i)
	{
		return 0;
	}
};




}
