// Created by Saliya Ekanayake on 1/22/19.


#ifndef LBL_DAL_ALPHABET_HPP
#define LBL_DAL_ALPHABET_HPP

#include <cstring>
#include <string>

#include "Types.hpp"

using std::string;	using std::cout;



class
Alphabet
{

public:

	enum type {PROTEIN, DNA};


	Alphabet (Alphabet::type t)
	{
		switch (t)
		{
		case PROTEIN:
			for (int i = 0; i < strlen(protein); ++i)
				al_map[protein[i]] = protein[i];
			break;

		case DNA:
			for (int i = 0; i < strlen(dna); ++i)
				al_map[dna[i]] = dna[i];
			break;
		}
	}

	

	void init (void);



	char &
	operator[] (size_t idx)
	{
		return letters[idx];
	}

	

	static const unsigned short capacity = 128;
	string						letters;
	unsigned short				size;
	unsigned short				max_char;
	unsigned char				char_to_code[capacity];
	unsigned char				code_to_char[capacity];
	unsigned char               al_map[capacity];


	static const char *protein;	// full superset
	static const char *dna;
};



class
DefaultProteinAlph : public Alphabet
{

public:

	DefaultProteinAlph () :
		Alphabet(Alphabet::PROTEIN)
	{
		letters = "ARNDCQEGHILKMFPSTWYVBZX*J";
		size	 = letters.size();
		max_char = 90;
		init();
	}		
};



class
Murphy10ProteinAlph : public Alphabet
{

public:

	Murphy10ProteinAlph () :
		Alphabet(Alphabet::PROTEIN)
	{
		letters = "ACDFGHIKSYBZX*J";

		al_map['D'] = al_map['E'] = al_map['N'] = al_map['Q'] = 'D';
		al_map['F'] = al_map['W'] = al_map['Y'] = 'F';
		al_map['I'] = al_map['L'] = al_map['M'] = al_map['V'] = 'I';
		al_map['K'] = al_map['R'] = 'K';
		al_map['S'] = al_map['T'] = 'S';		

		size	 = letters.size();
		max_char = 90;
		init();
	}		
};



class
DSSP10ProteinAlph : public Alphabet
{

public:

	DSSP10ProteinAlph () :
		Alphabet(Alphabet::PROTEIN)
	{
		letters = "EILFAWHCDGBZX*J";

		al_map['E'] = al_map['K'] = al_map['Q'] = al_map['R'] = 'E';
		al_map['I'] = al_map['V'] = 'I';		
		al_map['L'] = al_map['Y'] = 'L';
		al_map['A'] = al_map['M'] = 'A';
		al_map['H'] = al_map['T'] = 'H';
		al_map['D'] = al_map['N'] = al_map['S'] = 'D';
		al_map['G'] = al_map['P'] = 'G';
		
		size	 = letters.size();
		max_char = 90;
		init();
	}		
};



class
GBMR10ProteinAlph : public Alphabet
{

public:

	GBMR10ProteinAlph () :
		Alphabet(Alphabet::PROTEIN)
	{
		letters = "GDNAYHCTSPBZX*J";

		al_map['A'] = al_map['E'] = al_map['F'] = al_map['I'] = al_map['K'] =
			al_map['L'] = al_map['M'] = al_map['Q'] = al_map['R'] =
			al_map['V'] = al_map['W'] = 'A';
		
		size	 = letters.size();
		max_char = 90;
		init();
	}		
};



class
TD10ProteinAlph : public Alphabet
{

public:

	TD10ProteinAlph () :
		Alphabet(Alphabet::PROTEIN)
	{
		letters = "PGEDTHIWALBZX*J";

		al_map['E'] = al_map['K'] = al_map['R'] = al_map['Q'] = 'E';
		al_map['D'] = al_map['S'] = al_map['N'] = 'D';
		al_map['H'] = al_map['C'] = 'H';
		al_map['I'] = al_map['V'] = 'I';
		al_map['W'] = al_map['Y'] = al_map['F'] = 'W';
		al_map['L'] = al_map['M'] = 'L';
		
		size	 = letters.size();
		max_char = 90;
		init();
	}		
};



class
DiamondProteinAlph : public Alphabet
{

public:

    DiamondProteinAlph () :
		Alphabet(Alphabet::PROTEIN)
	{
		letters = "KCGHIMFYWPSBZX*J";

		al_map['K'] = al_map['R'] = al_map['E'] = al_map['D'] = al_map['Q'] =
			al_map['N'] = 'K';
		al_map['I'] = al_map['L'] = al_map['V'] = 'I';
		al_map['S'] = al_map['T'] = al_map['A'] = 'S';		
		
		size	 = letters.size();
		max_char = 90;
		init();
	}		
};



class
DefaultDNAAlph : public Alphabet
{

	DefaultDNAAlph () :
		Alphabet(Alphabet::DNA)
	{
		cout << "DefaultDNAAlph ctor\n";
		letters = "ACGT";
		size	 = letters.size();
		max_char = 84;
		init();
	}
};



#endif //LBL_DAL_ALPHABET_HPP
