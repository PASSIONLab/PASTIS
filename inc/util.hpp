/**
 * @file
 *	util.hpp
 *
 * @author
 *	Oguz Selvitopi
 *
 * @date
 *
 * @brief
 *	Simple util utility
 *
 * @todo
 *
 * @note
 *	
 */

#pragma once

#include <string>
#include <utility>
#include <vector>

#include "ParallelOps.hpp"
#include "Types.hpp"

extern shared_ptr<pastis::ParallelOps> parops;



namespace
pastis
{

std::string ptr_to_str (void *p);
std::string mb_str (uint64_t sz);
std::string gb_str (uint64_t sz);
std::string kb_str (uint64_t sz);



template <typename T>
std::vector<std::pair<T, T>>
cbmat_offsets
(
    T		sz,
	int		nblocks,
	bool	is_col
)
{
	int np	= parops->grid->GetGridRows();
	int pid = parops->grid->GetRankInProcCol();
	if (is_col)
	{
		np	= parops->grid->GetGridCols();
		pid = parops->grid->GetRankInProcRow();
	}

	T cur_range = 0;
	std::vector<std::pair<T, T>> out;
	for (int bi = 0; bi<nblocks; ++bi)
	{		
		T	bsz		= (sz + (nblocks-bi-1)) / nblocks;
		T	avg_len = bsz/np;
		T	len		= (pid == np-1) ? (bsz-(pid*avg_len)) : avg_len;
		T	offset	= pid*avg_len;
		out.emplace_back(cur_range+offset, len);
		cur_range += bsz;
	}

	return out;
}


 
}
