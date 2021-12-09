/**
 * @file
 *	sim_search.hpp
 *
 * @author
 *	Oguz Selvitopi
 *
 * @date
 *
 * @brief
 *	PASTIS pipeline for computing seq similarity matrix
 *
 * @todo
 *
 * @note
 *	
 */

#pragma once

#include "Types.hpp"



namespace
pastis
{



void compute_sim_mat_blocked (params_t &params);

template<typename T_OUT>
void
compute_sim_mat (params_t &params);

	
	
}
