//
// Created by Saliya Ekanayake on 2019-02-19.
// Modified by Aydin Buluc on 2019-12-29
//


#pragma once

#include <memory>

#include "DistFastaData.hpp"
#include "Types.hpp"



namespace
pastis
{


template <typename NT>
void
pw_aln_batch
(
 	std::shared_ptr<DistFastaData>		 dfd,
	combblas::SpParMat <uint64_t, NT, combblas::SpDCCols <uint64_t, NT>> &C,
	PWAlign 							&pwa,
	params_t							&params,
	int									 bri = -1,
	int									 bci = -1
);


template <typename NT>
void
pw_aln_batch_ovlp
(
 	std::shared_ptr<DistFastaData>		 dfd,
	combblas::SpParMat <uint64_t, NT, combblas::SpDCCols <uint64_t, NT>> &C,
	PWAlign 							&pwa,
	params_t							&params,
	MultData<MatrixEntry, NT>	    	&md,
	int									 bri = -1,
	int									 bci = -1
);

	
}
