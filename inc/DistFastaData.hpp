// Created by Saliya Ekanayake on 2019-02-17.

#pragma once

#include <string>
#include <utility>

#include "FastaData.hpp"
#include "ParallelOps.hpp"
#include "TraceUtils.hpp"
#include "Types.hpp"

extern shared_ptr<pastis::ParallelOps> parops;



namespace
pastis
{



// Stores information for sequences needed by a proc
struct
NbrData
{
	ushort		rc_flag;
	int			owner_rank;
	int			nbr_rank;
	int 		tag;			// message tag
	uint64_t	nbr_seq_start_idx;
	uint64_t	nbr_seq_end_idx;	

	NbrData () = default;

	NbrData (ushort rc_flag, int owner_rank, int nbr_rank, int tag,
			 uint64_t nbr_seq_start_idx, uint64_t nbr_seq_end_idx)
		: rc_flag(rc_flag),
		  owner_rank(owner_rank),
		  nbr_rank(nbr_rank),
		  tag(tag),
		  nbr_seq_start_idx(nbr_seq_start_idx),
		  nbr_seq_end_idx(nbr_seq_end_idx)
	{
	}



	std::string
	get_str (void)
	{
		std::string tmp = "row";
		if (rc_flag == 0)
			tmp = "col";
		return "neighbor rank " + std::to_string(nbr_rank) +
			" r/c " + tmp +
			" start,end " +
			(std::to_string(nbr_seq_start_idx) +
			 "," + std::to_string(nbr_seq_end_idx))
			+ " count " +
			std::to_string(nbr_seq_end_idx-nbr_seq_start_idx+1);
	}
};



class
DistFastaData
{

private:

	uint64_t				overlap;
	ushort					k;
	
	uint64_t	 l_seq_count;
  	uint64_t	*l_seq_counts	  = nullptr;
	uint64_t	 orig_g_seq_offset;
	uint64_t	 g_seq_count;
  	uint64_t	*g_seq_offsets	  = nullptr;
	bool		 is_diagonal_cell = false;	

	uint64_t row_seq_start_idx;
	uint64_t row_seq_end_idx;
	uint64_t col_seq_start_idx;
	uint64_t col_seq_end_idx;

	// @OGUZ-COMMENT move these to pairwise aligner
	// std::vector<seqan::Peptide *> row_seqs;
	// std::vector<seqan::Peptide *> col_seqs;

	// recv buffers
	int recv_nbrs_count;
	std::vector<int> recv_nbrs_idxs;
	char **recv_nbrs_buffs = nullptr;
	uint64_t *recv_nbrs_buff_lengths = nullptr;
	MPI_Request *recv_nbrs_buffs_reqs = nullptr;
	MPI_Status *recv_nbrs_buffs_stats = nullptr;
	// FastaData **recv_fds = nullptr;

	// send buffers
	int to_nbrs_count;
	MPI_Request *to_nbrs_buffs_reqs = nullptr;
	MPI_Status *to_nbrs_buffs_stat = nullptr;

	bool ready = false;			// seq comm over?


	// comm for blocked multiplication
	std::vector<std::pair<NbrData *, uint64_t>> bl_cnbrs_; // (nbr, recv size)	
	int bl_nsnbrs_;										   // number of sends
	MPI_Request *bl_rreqs_ = nullptr;
	MPI_Request *bl_sreqs_ = nullptr;
	vector<uint64_t> bl_rseq_offsets_;
	vector<uint64_t> bl_cseq_offsets_;
	

public:

	FastaData *fd = nullptr;
	std::vector<NbrData> my_nbrs;

	std::vector<NbrData> bl_nbrs_; // blocked multiplication
	std::vector<std::pair<char *, FastaData *>> bl_rbufs_; // (recv buf, fd)
		
	// received fasta data for grid seqs
	std::vector<FastaData *> recv_fds_;	



public:

	DistFastaData(const char *file, const std::string &idx_map_file,
				  uint64_t overlap, ushort k, params_t &params);

	~DistFastaData();



	uint64_t
	global_count ()
	{
		return g_seq_count;
	}



	uint64_t
	global_start_idx ()
	{
		return g_seq_offsets[parops->g_rank];
	}



	FastaData *
	lfd ()
	{
		return fd;
	}



	bool
	is_ready ()
	{
		return ready;
	}



	void
	wait ();



	void
	bl_wait ();




	uint64_t
	get_nrseqs ()
	{
		return row_seq_end_idx - row_seq_start_idx + 1;
	}



	uint64_t
	get_ncseqs ()
	{
		return col_seq_end_idx - col_seq_start_idx + 1;
	}



	uint64_t
	get_nrseqs_bl ()
	{
		return bl_rseq_offsets_[bl_rseq_offsets_.size()-1];
	}



	uint64_t
	get_ncseqs_bl ()
	{
		return bl_cseq_offsets_[bl_cseq_offsets_.size()-1];
	}



	uint64_t
	bl_rseq_offset (int b)
	{
		return bl_rseq_offsets_[b];
	}



	uint64_t
	bl_cseq_offset (int b)
	{
		return bl_cseq_offsets_[b];
	}




private:

	void
	read_fasta (const char *file, char *&buff, uint64_t &l_start,
				uint64_t &l_end);

	void
	write_idx_map (const char *idx_map_file);

	void
	collect_grid_seqs ();

	void
	collect_seqs (int br, int bc);

	void
	find_nbrs (int grid_rc_procs_count, int grid_rc_id,
			   uint64_t avg_l_seq_count, uint64_t avg_grid_seq_count,
			   uint64_t rc_seq_start_idx, ushort rc_flag);	
};

}
