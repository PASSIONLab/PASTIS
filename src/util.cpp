/**
 * @file
 *	util.cpp
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

#include <sstream>

#include "../inc/util.hpp"



namespace
pastis
{

std::string
ptr_to_str (void *p)
{
	std::ostringstream addr;
	addr << p;
	return addr.str();
}



std::string
mb_str (uint64_t sz)
{
	return std::to_string(sz/static_cast<double>(1<<20)) + " MB";
}



std::string
gb_str (uint64_t sz)
{
	return std::to_string(sz/static_cast<double>(1<<30)) + " GB";
}



std::string
kb_str (uint64_t sz)
{
	return std::to_string(sz/static_cast<double>(1<<10)) + " KB";
}



// template <class NT>
// uint64_t
// get_mat_size (PSpMat<NT>::MPI_DCCols &A)
// {
// 	return (A.getlocalnnz() * sizeof(uint64_t)) +
// 		(A.getlocalnnz() * sizeof(MatrixEntry)) +
// 		(A.getlocalcols() * sizeof(uint64_t) * 2);
// }

	
}
