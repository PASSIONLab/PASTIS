// Created by Saliya Ekanayake on 1/29/19.

#ifndef LBL_DAL_SPARSEMAT_HPP
#define LBL_DAL_SPARSEMAT_HPP

#include <cstdint>
#include <utility>
#include "CombBLAS/CombBLAS.h"

template <class NT>
class PSpMat
{
public:
  typedef combblas::SpTuples <uint64_t, NT> Tuples;	
  typedef combblas::SpDCCols <uint64_t, NT> DCCols;
  typedef combblas::SpParMat <uint64_t, NT, DCCols> MPI_DCCols;
  typedef std::tuple<uint64_t, uint64_t, NT *> ref_tuples;
};

#endif //LBL_DAL_SPARSEMAT_HPP
