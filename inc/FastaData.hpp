// Created by Saliya Ekanayake on 1/7/19.

#pragma once

#include <cstdlib>
#include <limits>
#include <memory>
#include <vector>

#include "macros.hpp"
#include "ParallelOps.hpp"
#include "TraceUtils.hpp"
#include "Types.hpp"

using std::shared_ptr;

extern shared_ptr<pastis::ParallelOps> parops;



namespace
pastis
{

class
FastaData
{

private:

	uvec_64		*id_starts	= nullptr;
  	uvec_64		*seq_starts = nullptr;
  	uvec_64		*del_idxs	= nullptr;
	char		*buff;
  	uint64_t	 l_seq_count;
	uint64_t	 l_start;
	uint64_t	 l_end;
	uint64_t	 orig_l_seq_count;
	



public:

	FastaData(char *buff, ushort k, uint64_t l_start, uint64_t &l_end);

	~FastaData();



	uint64_t local_count ()
	{
		return l_seq_count;
	}



	uint64_t orig_local_count ()
	{
		return orig_l_seq_count;
	}



	uvec_64 *
	deleted_indices ()
	{
		return del_idxs;
	}



	void
	buffer_size (uint64_t	 start_idx,
				 uint64_t	 end_idx_inclusive,
				 uint64_t	&len,
				 uint64_t	&start_offset,
				 uint64_t	&end_offset_inclusive)
	{
		assert(start_idx >= 0 && start_idx < id_starts->size() &&
			   end_idx_inclusive >= 0 && end_idx_inclusive < id_starts->size());
		start_offset = (*id_starts)[start_idx];
		end_offset_inclusive = (end_idx_inclusive+1 < id_starts->size() ?
								((*id_starts)[end_idx_inclusive+1] - 2) :
								l_end);
		len = end_offset_inclusive - start_offset + 1;
		#ifdef PASTIS_BOUND_CHECK
		assert(start_offset >= l_start &&
			   start_offset <= l_end &&
			   end_offset_inclusive >= l_start &&
			   end_offset_inclusive <= l_end);
		#endif
	}



	const char *
	buffer ()
	{
		return buff;
	}



	char *
	get_sequence (uint64_t idx, ushort &len, uint64_t &start_offset,
				  uint64_t &end_offset_inclusive)
	{
		// @OGUZ-WARNING trouble if sequence longer than max ushort
		uint64_t end_pos = idx+1 < id_starts->size() ?
			(*id_starts)[idx+1]-1 : l_end+1;

		assert(end_pos <= l_end+1);
		assert(idx < seq_starts->size());
		
		len = static_cast<ushort>(end_pos - (*seq_starts)[idx]);

		assert(len <= std::numeric_limits<ushort>::max());
		
		start_offset = (*seq_starts)[idx];
		end_offset_inclusive = end_pos-1;
		#ifdef PASTIS_BOUND_CHECK
		assert(start_offset >= l_start &&
			   start_offset <= l_end &&
			   end_offset_inclusive >= l_start &&
			   end_offset_inclusive <= l_end);
		#endif
		return buff;
	}



	uint64_t
	get_l_start ()
	{
		return l_start;
	}



	uint64_t
	get_l_end ()
	{
		return l_end;
	}
};

}
