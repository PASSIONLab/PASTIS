// Created by Saliya Ekanayake on 1/25/19.

#include <algorithm>
#include <cassert>
#include <cctype>
#include <iostream>

#include "../inc/Alphabet.hpp"

using std::fill_n;



namespace
pastis
{

/*!
 * NOTE: We have added the base 'J' as well.
 * Sometimes it is not possible two differentiate two closely related amino
 * acids, therefore we have the special cases:
 * asparagine/aspartic acid - asx - B
 * glutamine/glutamic acid - glx - Z
 * http://www.cryst.bbk.ac.uk/education/AminoAcid/the_twenty.html
 */
const char* Alphabet::protein = "ARNDCQEGHILKMFPSTWYVBZX*J";
const char* Alphabet::dna = "ACGT";



void
Alphabet::init (void)
{
	unsigned char no_code = capacity;
	fill_n(char_to_code, capacity, no_code);

	unsigned char code = 0;
	for (int i = 0; i < strlen(protein); ++i)
	{
		char c = protein[i];
		char d = al_map[c];
		// std::cout << c << ' ' << d << " | " << std::endl;
		if (char_to_code[d] == no_code)
		{
			char_to_code[d]	   = code;
			code_to_char[code] = static_cast<unsigned char>(d);
			++code;
		}
		char_to_code[c] = char_to_code[d];
	}
}

}
