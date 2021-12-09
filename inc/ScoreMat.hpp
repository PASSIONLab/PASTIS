// Created by Saliya Ekanayake on 2019-07-23.

#pragma once

#include <limits>
#include <map>
#include <string>
#include <vector>

#include "Types.hpp"



namespace
pastis
{

struct
Penalty
{
	
	char	base_char;
    ushort	sub_idx;
    short	penalty;


	
    // Penalty(char base_char);
    
};



struct
ScoreMatrix
{

	// max ascii char for protein alphabet is Z, which has the code 90
	ushort const row_size = 91;
    ushort const alph_size;
	char score_[91*91]{};
	std::map<char, std::vector<char>*> base_to_subtitutes; // subs for each char
	std::vector<Penalty> penalties[91]{};



	explicit
	ScoreMatrix (ushort alph_size) :
		alph_size(alph_size)
	{
	}



	inline
	short
	score (char ci, char cj)
	{
		return static_cast<short>(score_[ci*row_size+cj]);
	}



	inline
	short dist (char ci, char cj)
	{
		return static_cast<short>(score(ci, ci) - score(ci, cj));
	}



	short
	self_score (const std::string &kmer)
	{
		short sum = 0;
		for (const char &c : kmer)
			sum += score_[c*row_size+c];
		return sum;
	}



	short
	self_score (const char *ptr)
	{
		short sum = 0;
		for (; *ptr; ++ptr)
			sum += score_[(*ptr)*row_size + (*ptr)];
		return sum;
	}	
};



struct
Blosum62 : ScoreMatrix
{

private:

	// chars only in protein alphabet (* and A-Z)
	char const data[25*25] =
		{
         4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-2,-1, 0,-4,-1,
        -1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-1, 0,-1,-4,-2,
        -2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 3, 0,-1,-4,-3,
        -2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, 4, 1,-1,-4,-4,
         0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4,-1,
        -1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0, 3,-1,-4,-2,
        -1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4,-3,
         0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-1,-2,-1,-4,-4,
        -2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, 0, 0,-1,-4,-3,
        -1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-3,-3,-1,-4, 2,
        -1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4,-3,-1,-4, 4,
        -1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, 0, 1,-1,-4,-2,
        -1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-3,-1,-1,-4, 2,
        -2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-3,-3,-1,-4, 0,
        -1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2,-1,-2,-4,-3,
         1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 0, 0, 0,-4,-2,
         0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,-1,-1, 0,-4,-1,
        -3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-3,-2,-4,-2,
        -2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-3,-2,-1,-4,-1,
         0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-3,-2,-1,-4, 1,
        -2,-1, 3, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4, 1,-1,-4,-4,
        -1, 0, 0, 1,-3, 3, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4,-3,
         0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-1,-1,-4,-1,
        -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1,-4,
        -1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4,-3,-1,-4, 4
		};



public:

	Blosum62 () :
		ScoreMatrix(25)
	{
		const char *alph = Alphabet::protein;
		for (int i = 0; i < alph_size; ++i)
		{
			char	ci			 = alph[i];
			ushort	score_offset = ci * row_size;
			ushort	data_offset	 = i * alph_size;

			std::vector<char> *subs = new std::vector<char>(alph_size);
			short max_score = std::numeric_limits<short>::min();
			short max_score_count = 0;

			// find max score among other chars and record them in subs
			for (int j = 0; j < alph_size; ++j)
			{
				char cj = alph[j];
				char s = data[data_offset + j];

				score_[score_offset + cj] = s;

				if (ci != cj)
				{
					if (s > max_score)
					{
						max_score_count = 0;
						max_score = s;
						(*subs)[max_score_count] = cj;
						++max_score_count;
					}
					else if (s == max_score)
					{
						(*subs)[max_score_count] = cj;
						++max_score_count;
					}
				}
			}

			(*subs)[max_score_count] = ci; // myself
			++max_score_count;
			subs->erase(subs->begin()+max_score_count, subs->end());
			base_to_subtitutes[ci] = subs;
		}
	}
};

}
