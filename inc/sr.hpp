/**
 * @file
 *  sr.hpp
 *
 * @author
 *  Oguz Selvitopi
 *
 * @date
 *
 * @brief
 *  Semirings used in matrix operations
 *
 * @todo
 *
 * @note
 *
 */

#pragma once

#include <algorithm>



namespace
pastis
{

// semiring for detecting overlap - kmer location version
template<typename IN,
		 typename OUT>
struct
KmerIntersectLoc
{

	static
	OUT
	id()
	{
		return OUT();
	}



	static
	bool
	returnedSAID ()
	{
		return false;
	}



	static
	OUT
	multiply (const IN &arg1, const IN &arg2)
	{
		OUT a;
		a.first.first  = arg1.offset;
		a.first.second = arg2.offset;
		a.score		   = std::max(arg1.cost, arg2.cost);
		return a;
    }



    static
	OUT
	add (const OUT &arg1, const OUT &arg2)
	{
		OUT res(arg1.count + arg2.count);		
		res.first.first	  = arg1.first.first;
		res.first.second  = arg1.first.second;
		res.second.first  = arg2.first.first;
		res.second.second = arg2.first.second;
		res.score		  = std::max(arg1.score, arg2.score);
		return res;
    }



	static
	void
	axpy (IN a, const IN &x, OUT &y)
	{
      	y = add(y, multiply(a, x));
    }



	static
	MPI_Op
	mpi_op()
	{
		static MPI_Op mpiop;
		static bool exists = false;
		if (exists)
			return mpiop;
		else
		{
			MPI_Op_create(MPI_func, true, &mpiop);
			exists = true;
			return mpiop;
		}
    }



	// @OGUZ-FIXED  *((OUT) inoutvec + 1) should be  *((OUT) inoutvec + i)
	static
	void
    MPI_func (void *invec, void *inoutvec, int *len, MPI_Datatype *datatype)
	{
      	for (int i = 0; i < *len; ++i)
			*((OUT)inoutvec+i) = add(*((OUT)invec+i), *((OUT)inoutvec+i));
    }
};



// semiring for detecting overlap - light version
template<typename IN,
		 typename OUT>
struct
KmerIntersectLight
{

	static
	OUT
	id()
	{
		return OUT();
	}



	static
	bool
	returnedSAID ()
	{
		return false;
	}



	static
	OUT
	multiply (const IN &arg1, const IN &arg2)
	{
		OUT a;
		a.score = std::max(arg1.cost, arg2.cost);
		return a;
    }



    static
	OUT
	add (const OUT &arg1, const OUT &arg2)
	{
		OUT res(arg1.count + arg2.count);
		res.score = std::max(arg1.score, arg2.score);
		return res;
    }



	static
	void
	axpy (IN a, const IN &x, OUT &y)
	{
      	y = add(y, multiply(a, x));
    }



	static
	MPI_Op
	mpi_op()
	{
		static MPI_Op mpiop;
		static bool exists = false;
		if (exists)
			return mpiop;
		else
		{
			MPI_Op_create(MPI_func, true, &mpiop);
			exists = true;
			return mpiop;
		}
    }



	// @OGUZ-FIXED  *((OUT) inoutvec + 1) should be  *((OUT) inoutvec + i)
	static
	void
    MPI_func (void *invec, void *inoutvec, int *len, MPI_Datatype *datatype)
	{
      	for (int i = 0; i < *len; ++i)
			*((OUT)inoutvec+i) = add(*((OUT)invec+i), *((OUT)inoutvec+i));
    }
};

}
