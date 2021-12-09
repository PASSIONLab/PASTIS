// Created by Saliya Ekanayake on 10/15/19.

#pragma once

#include <unordered_set>
#include <vector>

#include "../Alphabet.hpp"
#include "../DistFastaData.hpp"
#include "../ParallelOps.hpp"
#include "../Types.hpp"
#include "Kmer.hpp"


namespace
pastis
{

PSpMat<MatrixEntry>::MPI_DCCols
generate_A (uint64_t seq_count, std::shared_ptr<DistFastaData> &dfd, ushort k,
			ushort s, Alphabet &alph,
			std::unordered_set<Kmer, Kmer> &local_kmers);

}
