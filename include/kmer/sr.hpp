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



namespace
pastis
{

template <typename IN,
		  typename OUT>
struct
CommonKmerCnt
{
	static
	OUT
	id ()
	{
		return 0;
	}



	static
	bool
	returnedSAID ()
	{
		return false;
	}



	static
	OUT
	add (const OUT &arg1, const OUT &arg2)
	{
		return arg1 + arg2;
	}



	static
	OUT
	multiply (const IN &arg1, const IN &arg2)
	{
		return 1;
	}



	static
	void
	axpy (IN a, const IN &x, OUT &y)
	{
		y = add(y, multiply(a, x));
	}
};
}
