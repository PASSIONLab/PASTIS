// Created by Saliya Ekanayake on 10/17/19.

#pragma once

#include <set>
#include <string>
#include <utility>
#include <vector>

#include "../Alphabet.hpp"



namespace
pastis
{

struct
Kmer
{

private:

	uint64_t	kmer_code;
	std::string	kmer_str;
	std::set<ushort, std::greater<ushort>> free_idxs;
	short		dist2r = 0;



public:

	Kmer ()
	{
	}



	Kmer (uint64_t kmer_code, ushort k, Alphabet &alph) :
		kmer_code(kmer_code)
	{
		kmer_str = std::move(std::string(k, ' '));
		uint64_t q, r;
		ushort free_idx = 0;
		ushort i = k-1;
		while (i >= 1)
		{
			q = kmer_code / alph.size;
			r = kmer_code - (q*alph.size);
			kmer_str[i] = alph.code_to_char[r];
			free_idxs.insert(free_idx++);
			kmer_code = q;
			--i;
		}
		kmer_str[i] = alph.code_to_char[kmer_code];
		free_idxs.insert(free_idx);
	}



	// @OGUZ-COMMENT currently not used
	// Kmer (std::string str, uint64_t kmer_code, Alphabet &alph, bool fix_kcode) :
	// 	kmer_code(kmer_code), kmer_str(std::move(str))
	// {
	// 	ushort free_idx = 0;
	// 	size_t len = kmer_str.length();
	// 	for (size_t i = 0; i < len; ++i)
	// 		free_idxs.insert(free_idx++);
	// 	if (fix_kcode)
	// 		update_kmer_code(alph);
	// }


	
	~Kmer ()
	{
	}



	// for use as hash function
	size_t
	operator() (const Kmer &t) const
    {
      	return t.kmer_code;
    }



	bool
	operator== (const Kmer& t) const
    {
      	return kmer_code == t.kmer_code;
    }



	inline
	uint64_t code () const
	{
		return kmer_code;
	}




private:

	void
	update_kmer_code (Alphabet &alph)
	{
		ushort base = alph.size;
		kmer_code = 0;
		for (char cap_c : kmer_str)
		{
			if (cap_c > 96 && cap_c < 123) // lower to upper case
				cap_c = cap_c - 32;
			kmer_code = kmer_code*base + alph.char_to_code[cap_c];
		}
    }
};

}
