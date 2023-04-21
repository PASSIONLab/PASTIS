// Created by Saliya Ekanayake on 1/7/19.

#include "../inc/FastaData.hpp"
#include "../inc/macros.hpp"

using std::to_string;


namespace
pastis
{

FastaData::FastaData
(
    char							*buff,
	ushort							 k,
	uint64_t						 l_start,
	uint64_t						&l_end
)
{
	id_starts	= new uvec_64();
  	seq_starts	= new uvec_64();
  	del_idxs	= new uvec_64();
	l_seq_count = 0;

	char	c;
  	bool	in_name = false;
  	bool	in_seq	= false;
  	int		seq_len = 0;
	// nc_count means No character count. This includes new line and *
	// characters It also includes entire sequences that are less than k-mer
	// length
  	uint64_t nc_count = 0;
  	uint64_t idx;

	// updates buff by overwriting deleted stuff. idx comes from behind by
	// writing valid characters onto buff
	for (uint64_t i = l_start; i <= l_end; ++i)
	{
		c		  = buff[i];
		idx		  = i - nc_count;
		buff[idx] = c;

		// Modify 'u' or 'U' in sequences to 'T'.
    	// This is according to http://meme-suite.org/doc/alphabets.html
		if ((c == 'U' || c == 'u') && in_seq)
		{
			buff[idx] = 'T';
		}

		// !in_name logic is important as some fasta files have > character in
		// the middle of the sequence identifier. If we didn't have this test
		// then we'll consider that as the start of another sequence
		if (c == '>' && !in_name) // start a new sequence identifier
		{
			id_starts->push_back(idx);
      		seq_len = 0;
      		in_name = true;
      		in_seq	= false;
      		++l_seq_count;
		}
		else if (c == '\n')		// end a seq identifier or seq
		{
			if (in_name && i+1 <= l_end) // end seq identifier
			{
				seq_starts->push_back(idx+1);
				in_name = false;
				in_seq	= true;
			}
			else if (in_seq && i+1 <= l_end) // end seq
			{
				if (buff[i+1] != '>')
				{
					++nc_count;
				}
				else if (buff[i+1] == '>')
				{
					if (seq_len < k) // remove sequence if smaller than kmer
									 // length
					{
						++seq_len; // capture the new line character too for removal
						uint64_t seq_id_start = id_starts->back();
						uint64_t seq_start = seq_starts->back();
            			// + 1 to capture the new line between sequence id and
            			// sequence start
						nc_count += seq_len + (seq_start-seq_id_start) + 1;
						id_starts->pop_back();
						seq_starts->pop_back();
						del_idxs->push_back(l_seq_count + del_idxs->size());
						--l_seq_count;
					}
				}
			}
		}
		else if (c == '*')
		{
			if (in_seq)
			{
				++nc_count;
			}
		} else
		{
			if (in_seq)
			{
				++seq_len;
			}
		}
	}


    // Remove the last sequence as well if it's shorter than k
	if (seq_len < k)
	{
    	// The last sequence doesn't have a new line at the end as our l_end
    	// points to the last real character of the block.
    	uint64_t seq_id_start = id_starts->back();
    	uint64_t seq_start = seq_starts->back();
    	// + 1 to capture the new line between sequence id and sequence start
    	nc_count += seq_len + (seq_start-seq_id_start) + 1;
    	id_starts->pop_back();
    	seq_starts->pop_back();
    	--l_seq_count;
	}


	assert(l_seq_count > 0);

	uint64_t tmp = l_end;

	this->buff		  = buff;
	this->l_start	  = l_start;
	l_end			 -= nc_count;
	this->l_end		  = l_end;
	orig_l_seq_count  = l_seq_count + del_idxs->size();

	parops->logger->log(string("local seq count ") + to_string(l_seq_count) +
						string(" eliminated ") + to_string(del_idxs->size()));
	parops->logger->log(string("l_start ") + to_string(l_start) +
						string(" l_end ") + to_string(l_end) +
						string(" nc_count ") + to_string(nc_count) +
						string(" orig l_end ") + to_string(tmp));

	#if PASTIS_DBG_LVL > 0
	parops->bytes_alloc += id_starts->capacity() * sizeof(uint64_t);
	parops->bytes_alloc += seq_starts->capacity() * sizeof(uint64_t);
	parops->bytes_alloc += del_idxs->capacity() * sizeof(uint64_t);
	#endif
}



FastaData::~FastaData ()
{
  	delete id_starts;
	delete seq_starts;
	delete del_idxs;
	if (buff != NULL)
	{
		// free(buff);
		delete buff;
	}
}

}
