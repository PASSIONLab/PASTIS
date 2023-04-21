// Created by Saliya Ekanayake on 2019-02-17.

#include <algorithm>
#include <cassert>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "../inc/DistFastaData.hpp"
#include "../inc/macros.hpp"
#include "../inc/util.hpp"

using std::string;	using std::vector;	using std::sort;
using std::upper_bound;	using std::for_each;	using std::pair;
using std::to_string;  




namespace
pastis
{

DistFastaData::DistFastaData
(
    const char						*file,
	const string					&idx_map_file,
	uint64_t						 overlap,
	ushort							 k,
	params_t						&params
) :
	overlap(overlap), k(k)
{
	is_diagonal_cell =
		parops->grid->GetRankInProcRow() == parops->grid->GetRankInProcCol();

	parops->tp->start_timer("fasta|io-seqs");
	
	char *buff = nullptr;
	uint64_t l_start, l_end;
	read_fasta(file, buff, l_start, l_end);
	
	parops->tp->stop_timer("fasta|io-seqs");

	parops->tp->start_timer("fasta|local-setup");
	
	fd = new FastaData(buff, k, l_start, l_end);
	l_seq_count = fd->local_count();
	
	parops->tp->stop_timer("fasta|local-setup");

	l_seq_counts = new uint64_t[parops->g_np];
	MPI_Allgather(&l_seq_count, 1, MPI_UINT64_T, l_seq_counts,
				  1, MPI_UINT64_T, MPI_COMM_WORLD);

	g_seq_offsets = new uint64_t[parops->g_np+1];	
	g_seq_offsets[0] = 0;
	for (int i = 1; i < parops->g_np; ++i)
		g_seq_offsets[i] = g_seq_offsets[i-1] + l_seq_counts[i-1];

	g_seq_count = g_seq_offsets[parops->g_np-1] + l_seq_counts[parops->g_np-1];
	g_seq_offsets[parops->g_np] = g_seq_count;

	uint64_t orig_l_seq_count = fd->orig_local_count();
	orig_g_seq_offset = 0;
	MPI_Exscan(&orig_l_seq_count, &orig_g_seq_offset, 1,
			   MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

	if (!idx_map_file.empty())
	{
		parops->logger->log("Writing sequence indices to " + idx_map_file);
		this->write_idx_map(idx_map_file.c_str());
	}

	string s_tmp;
	parops->logger->log("l_seq_count " + to_string(l_seq_count));
	parops->logger->log("orig_l_seq_count " + to_string(orig_l_seq_count));	
	s_tmp.append("l_seq_counts ");
	for (int i = 0; i < parops->g_np; ++i)
		s_tmp.append(to_string(l_seq_counts[i]) + " ");
	parops->logger->log(s_tmp);
	s_tmp = "";
	s_tmp.append("g_seq_offsets ");
	for (int i = 0; i < parops->g_np+1; ++i)
		s_tmp.append(to_string(g_seq_offsets[i]) + " ");
	parops->logger->log(s_tmp);

	parops->logger->log("g_seq_count " + to_string(g_seq_count));
	parops->logger->log("orig_g_seq_offset " + to_string(orig_g_seq_offset));

	parops->tp->start_timer("fasta|seq-comm-setup");

	if (params.br == -1)
		collect_grid_seqs();
	else
		collect_seqs(params.br, params.bc);

	parops->tp->stop_timer("fasta|seq-comm-setup");
}



DistFastaData::~DistFastaData ()
{
	// for (auto &row_seq : row_seqs)
	// 	delete (row_seq);

  	// if (!is_diagonal_cell)
	// {
	//   	/*! If this was a diagonal cell then both the row_seqs and col_seqs
	// 	 * would point to the same sequences, so deleting row_seqs is enough.
	// 	 */
	// 	for (auto &col_seq : col_seqs)
	// 		delete (col_seq);
	// }

	// if (recv_fds != nullptr)
	// {
	// 	for (int i = 0; i < recv_nbrs_count; ++i)
	// 		delete (recv_fds[i]);
	// 	delete []recv_fds;
	// }
	
	delete[] to_nbrs_buffs_stat;
	delete[] to_nbrs_buffs_reqs;
	delete[] recv_nbrs_buffs_stats;
	delete[] recv_nbrs_buffs_reqs;
	
	delete[] recv_nbrs_buff_lengths;
	// // delete (fd->buffer());
	delete fd;
	delete[] l_seq_counts;
	delete[] g_seq_offsets;

	for (FastaData *fd : recv_fds_)
		delete fd;

	if (recv_nbrs_buffs != nullptr)
	{
		// below freed by FastaData
		// for (int i = 0; i < recv_nbrs_count; ++i)
		// {
		// 	if (recv_nbrs_buffs[i] != nullptr)
		// 		delete[] recv_nbrs_buffs[i];
		// }
		delete[] recv_nbrs_buffs;
	}

	if (bl_rreqs_)
		delete[] bl_rreqs_;
	if (bl_sreqs_)
		delete[] bl_sreqs_;
	for (auto &p : bl_rbufs_)
		delete p.second;		// p.first is freed by FastaData
		
}



void
DistFastaData::read_fasta
(
    const char	 *file,
	char		*&buff,
	uint64_t	 &l_start,
	uint64_t	 &l_end
)
{
	MPI_File f;
	int err = MPI_File_open(MPI_COMM_WORLD, file,
							MPI_MODE_RDONLY, MPI_INFO_NULL, &f);
	if (err)
	{
		if (parops->g_rank == 0)
			fprintf(stderr, "Couldn't open file %s\n", file);
		parops->teardown_parallelism();
		exit(1);
  	}

	int rank = parops->g_rank, np = parops->g_np;
	MPI_Offset g_start, g_end, file_size;
  	uint64_t l_chunk_size;
	double tmp;

	/* figure out who reads what */
  	MPI_File_get_size(f, &file_size);
  	file_size--;  /* get rid of text file eof */
  	l_chunk_size = static_cast<uint64_t >(file_size/np);
	g_start = rank * l_chunk_size;
	g_end = g_start + l_chunk_size - 1;
  	if (rank == np-1)
		g_end = file_size-1;

	/* add overlap to the end of everyone's chunk except last proc */
  	if (rank != np-1)
		g_end += overlap;

	buff = nullptr;
	l_chunk_size = static_cast<uint64_t >(g_end - g_start + 1);
	// buff = static_cast<char *>(malloc((l_chunk_size+1) * sizeof(char)));
	// buff = static_cast<char *>(malloc((l_chunk_size+1) * sizeof(char)));
	buff = new char[l_chunk_size+1];
	assert(buff != nullptr);

	#if PASTIS_DBG_LVL > 0
	parops->bytes_alloc += (l_chunk_size+1) * sizeof(char);
	#endif
		
	tmp = file_size;
	parops->logger->log(string("Parallel reading fasta file, size = ") +
						to_string(tmp) + string(" B, ") +
						to_string(tmp/(1<<10)) + string(" KB, ") +
						to_string(tmp/(1<<20)) + string(" MB, ") +
						to_string(tmp/(1<<30)) + string(" GB"));
	parops->logger->log(string("chunk size ") + to_string(l_chunk_size) +
						string(", read beg " ) + to_string(g_start) +
						string(", read end " ) + to_string(g_end));
	tmp = (l_chunk_size+1) * sizeof(char);
	parops->logger->log(string("buff size ") +
						to_string(tmp/(1<<10)) + string(" KB, ") +
						to_string(tmp/(1<<20)) + string(" MB, ") +
						to_string(tmp/(1<<30)) + string(" GB"));
	// string s_tmp = "Addresses: "
	// 	"f (MPI_File) " + ptr_to_str(static_cast<void *>(f)) +
	// 	" buff beg (char *) " + ptr_to_str(buff) +
	// 	" buff end (char *) " + ptr_to_str(buff+l_chunk_size);		
	// parops->logger->log(s_tmp);

	// @OGUZ-TODO read multiple times if l_chunk_size > max int
	assert(l_chunk_size <= std::numeric_limits<int>::max());
	MPI_Status status;
	MPI_File_read_at_all(f, g_start, buff, static_cast<int>(l_chunk_size),
						 MPI_CHAR, &status);
	buff[l_chunk_size] = '\0';

	int count;
	MPI_Get_count(&status, MPI_CHAR, &count);
	assert(count == l_chunk_size);
	parops->logger->log("MPI_File_read_at_all status: count " +
						to_string(count));

	#if PASTIS_DBG_LVL > 1
	parops->logger->log("ensuring contents of buff");
	for (uint64_t i = 0; i < l_chunk_size; ++i)
		assert(buff[i] > 31 && buff[i] < 127);
	#endif

	l_start = 0, l_end = l_chunk_size - 1;
  	if (rank != 0)
	{
		while (buff[l_start] != '>')
			++l_start;
	}

	if (rank != np-1)
	{
		l_end -= overlap;
		while (buff[l_end] != '>')
			++l_end;
		l_end -= 2;					// don't need > and '\n' at the end
  	}
}



void
DistFastaData::write_idx_map
(
	const char *idx_map_file 
)
{
	uint64_t			 g_seq_offset = g_seq_offsets[parops->g_rank];
	uvec_64				*del_idxs	  = fd->deleted_indices();
	uint64_t			 del_count	  = del_idxs->size();
	uint64_t			 cur_dels	  = 0;
    uint64_t			 orig_idx;
	std::stringstream	 ss;

	for (uint64_t i = 0; i < l_seq_count; ++i)
	{
		// this index was deleted
		if ((cur_dels<del_count) && (i+cur_dels == (*del_idxs)[cur_dels]))
			++cur_dels;

		orig_idx = i + cur_dels + orig_g_seq_offset;
		ss << (i+g_seq_offset) << " " << orig_idx << "\n";
	}

	parops->write_file_in_parallel(idx_map_file, ss.str());
}



void
DistFastaData::collect_grid_seqs ()
{
	int			pr				   = parops->grid->GetGridRows();
  	int			pc				   = parops->grid->GetGridCols();
	int			grid_row_id		   = parops->grid->GetRankInProcCol();
  	int			grid_col_id		   = parops->grid->GetRankInProcRow();
	uint64_t	avg_l_seq_count	   = g_seq_count / parops->g_np;
  	uint64_t	avg_grid_seq_count = g_seq_count / pr;

	row_seq_start_idx = avg_grid_seq_count * grid_row_id;
  	row_seq_end_idx = (grid_row_id == pr - 1) ? g_seq_count :
		row_seq_start_idx + avg_grid_seq_count;
	--row_seq_end_idx;			// get the zero based inclusive end index

	col_seq_start_idx = avg_grid_seq_count * grid_col_id;
  	col_seq_end_idx = (grid_col_id == pc - 1) ? g_seq_count :
		col_seq_start_idx + avg_grid_seq_count;
  	--col_seq_end_idx;			// get the zero based inclusive end index


	string s_tmp = "grind rank " + to_string(grid_row_id) +
		", " + to_string(grid_col_id) +
		" avg local seq count: " + to_string(avg_l_seq_count) +
		" avg grid seq count: " + to_string(avg_grid_seq_count) +
		" row seq range (grid): [" + to_string(row_seq_start_idx) +
		", " + to_string(row_seq_end_idx) + "]" +
		" col seq range (grid): [" + to_string(col_seq_start_idx) +
		", " + to_string(col_seq_end_idx) + "]";
	parops->logger->log(s_tmp);

	parops->logger->log("finding neighbors for seqs");
	

	find_nbrs(pr, grid_row_id, avg_l_seq_count, avg_grid_seq_count,
			  row_seq_start_idx, 1);
	if (grid_row_id != grid_col_id)
		find_nbrs(pc, grid_col_id, avg_l_seq_count, avg_grid_seq_count,
				  col_seq_start_idx, 0);

	this->recv_nbrs_count = 0;
	for (int i = 0; i < my_nbrs.size(); ++i)
	{
		if (my_nbrs[i].nbr_rank == parops->g_rank)
			continue;
		this->recv_nbrs_idxs.push_back(i);
		++(this->recv_nbrs_count);
	}

	parops->logger->log("getting buffer sizes for seqs");

	// get sizes for the buffers of the seqs I'll receive
	auto *recv_nbrs_buff_lengths_reqs = new MPI_Request[this->recv_nbrs_count];
	this->recv_nbrs_buff_lengths = new uint64_t[this->recv_nbrs_count];
	for (int i = 0; i < this->recv_nbrs_count; ++i)
		MPI_Irecv(this->recv_nbrs_buff_lengths + i,
				  1,
				  MPI_UINT64_T,
				  my_nbrs[recv_nbrs_idxs[i]].nbr_rank,
				  99+my_nbrs[recv_nbrs_idxs[i]].rc_flag,
				  MPI_COMM_WORLD,
				  recv_nbrs_buff_lengths_reqs + i);


	// who will need my seqs
	int block_length[6] = {1, 1, 1, 1, 1, 1};
	MPI_Aint displacement[6] = {offsetof(NbrData, rc_flag),
								offsetof(NbrData, owner_rank),
								offsetof(NbrData, nbr_rank),
								offsetof(NbrData, tag),
								offsetof(NbrData, nbr_seq_start_idx),
								offsetof(NbrData, nbr_seq_end_idx)};
	MPI_Datatype types[] = {MPI_UNSIGNED_SHORT, MPI_INT, MPI_INT, MPI_INT,
							MPI_UINT64_T, MPI_UINT64_T};
	MPI_Datatype MPI_NbrData;
	MPI_Type_create_struct(6, block_length, displacement, types, &MPI_NbrData);
	MPI_Type_commit(&MPI_NbrData);

	parops->logger->log("communicating neighbor metadata");

	// gather all neighbor metadata first
	int my_nbrs_count = static_cast<int>(my_nbrs.size());
	auto *all_nbrs_counts = new int[parops->g_np];
	MPI_Allgather(&my_nbrs_count,
				  1,
				  MPI_INT,
				  all_nbrs_counts,
				  1,
				  MPI_INT,
				  MPI_COMM_WORLD);

	int all_nbrs_count = 0;
	auto *all_nbrs_displas = new int[parops->g_np];
	all_nbrs_displas[0] = 0;
	for (int i = 0; i < parops->g_np; ++i)
	{
		all_nbrs_count += all_nbrs_counts[i];
		if (i > 0)
			all_nbrs_displas[i] = all_nbrs_displas[i-1] + all_nbrs_counts[i-1];
	}

	auto *all_nbrs = new NbrData[all_nbrs_count];
	MPI_Allgatherv(&my_nbrs[0],
				   my_nbrs_count,
				   MPI_NbrData,
				   all_nbrs,
				   all_nbrs_counts,
				   all_nbrs_displas,
				   MPI_NbrData,
				   MPI_COMM_WORLD);

	vector<uint64_t> send_lengths, send_start_offsets;
	vector<int> to_nbrs_idxs;
	to_nbrs_count = 0;
	sort(all_nbrs, all_nbrs + all_nbrs_count,
		 [](const NbrData &a, const NbrData &b) -> bool
		 {
			 return a.nbr_rank < b.nbr_rank;
		 });
	s_tmp = "";
	for (int i = 0; i < all_nbrs_count; ++i)
	{
		if (all_nbrs[i].nbr_rank != parops->g_rank) // doesn't request from me
			continue;
		if (all_nbrs[i].owner_rank == parops->g_rank) // local seqs
			continue;

		NbrData &nbr = all_nbrs[i];
		uint64_t len, start_offset, end_offset;
		fd->buffer_size(nbr.nbr_seq_start_idx, nbr.nbr_seq_end_idx,
						len, start_offset, end_offset);
		send_lengths.push_back(len);
		send_start_offsets.push_back(start_offset);
		to_nbrs_idxs.push_back(i);
		++to_nbrs_count;

		s_tmp += "to P " + to_string(all_nbrs[i].owner_rank) +
			" " + to_string(len/static_cast<double>(1<<20)) + " | ";
	}

	parops->logger->log("Seq send sizes (MB) " + s_tmp);

	// send lengths of the seqs to be transferred - matches Irecv's above
	auto *to_nbrs_send_reqs = new MPI_Request[to_nbrs_count];
	for (int i = 0; i < to_nbrs_count; ++i)
		MPI_Isend(&send_lengths[i],
				  1,
				  MPI_UINT64_T,
				  all_nbrs[to_nbrs_idxs[i]].owner_rank,
				  99+all_nbrs[to_nbrs_idxs[i]].rc_flag,
				  MPI_COMM_WORLD,
				  to_nbrs_send_reqs+i);

	// Wait for length transfers
	auto recv_stats = new MPI_Status[recv_nbrs_count];
	auto send_stats = new MPI_Status[to_nbrs_count];
	MPI_Waitall(recv_nbrs_count, recv_nbrs_buff_lengths_reqs, recv_stats);
	MPI_Waitall(to_nbrs_count, to_nbrs_send_reqs, send_stats);

	s_tmp = "";
	for (int i = 0; i < recv_nbrs_count; ++i)
		s_tmp += "from P " + to_string(my_nbrs[recv_nbrs_idxs[i]].nbr_rank) +
			" " + to_string(recv_nbrs_buff_lengths[i]/
						  static_cast<double>(1<<20)) + " | ";
	parops->logger->log("Seq recv sizes (MB) " + s_tmp);

	
	for (int i = 0; i < to_nbrs_count; ++i)
	{
		if (send_lengths[i] > std::numeric_limits<int>::max())
			parops->logger->log("Send seq size larger than INT_MAX: " +
								to_string(send_lengths[i]) + " bytes from P " +
								to_string(parops->g_rank) + " to P " +
								to_string(all_nbrs[to_nbrs_idxs[i]].owner_rank),
								Logger::LogLevel::WARNING);
	}
	for (int i = 0; i < recv_nbrs_count; ++i)
	{
		if (recv_nbrs_buff_lengths[i] > std::numeric_limits<int>::max())
			parops->logger->log("Recv seq size larger than INT_MAX: " +
								to_string(recv_nbrs_buff_lengths[i]) +
								" bytes from P " +
								to_string(my_nbrs[recv_nbrs_idxs[i]].nbr_rank) +
								" to P " +
								to_string(parops->g_rank),
								Logger::LogLevel::WARNING);
	}		

	// @OGUZ-TODO do multiple send/recvs below if buff length > int max


	parops->logger->log("starting communication of seqs");

	#if PASTIS_DBG_LVL > 0
	for (int i = 0; i < recv_nbrs_count; ++i)
		parops->bytes_alloc += recv_nbrs_buff_lengths[i] * sizeof(char);
	#endif
	
	// communicate seqs
	recv_nbrs_buffs = new char *[recv_nbrs_count];
	for (int i = 0; i < recv_nbrs_count; ++i)
		recv_nbrs_buffs[i] = new char[recv_nbrs_buff_lengths[i]];
	recv_nbrs_buffs_reqs = new MPI_Request[recv_nbrs_count];
  	recv_nbrs_buffs_stats = new MPI_Status[recv_nbrs_count];
	for (int i = 0; i < recv_nbrs_count; ++i)
		MPI_Irecv(recv_nbrs_buffs[i],
				  static_cast<int>(recv_nbrs_buff_lengths[i]),
				  MPI_CHAR,
				  my_nbrs[recv_nbrs_idxs[i]].nbr_rank,
				  77+my_nbrs[recv_nbrs_idxs[i]].rc_flag,
				  MPI_COMM_WORLD,
				  recv_nbrs_buffs_reqs+i);

	to_nbrs_buffs_reqs = new MPI_Request[to_nbrs_count];
  	to_nbrs_buffs_stat = new MPI_Status[to_nbrs_count];
	for (int i = 0; i < to_nbrs_count; ++i)
		MPI_Isend(fd->buffer() + send_start_offsets[i],
				  static_cast<int>(send_lengths[i]),
				  MPI_CHAR,
				  all_nbrs[to_nbrs_idxs[i]].owner_rank,
				  77+all_nbrs[to_nbrs_idxs[i]].rc_flag,
				  MPI_COMM_WORLD,
				  to_nbrs_buffs_reqs+i);


	delete[] recv_nbrs_buff_lengths_reqs;
	delete[] all_nbrs_counts;
	delete[] all_nbrs_displas;
	delete[] all_nbrs;
	delete[] to_nbrs_send_reqs;
	MPI_Type_free(&MPI_NbrData);
	delete[] recv_stats;
	delete[] send_stats;
}



void
DistFastaData::find_nbrs
(
    int			grid_rc_procs_count,
	int			grid_rc_id,
	uint64_t	avg_l_seq_count,
	uint64_t	avg_grid_seq_count,
	uint64_t	rc_seq_start_idx,
	ushort		rc_flag
)
{
	int start_rank = static_cast<int>((rc_seq_start_idx+1)/avg_l_seq_count);
	while (rc_seq_start_idx < g_seq_offsets[start_rank])
		--start_rank;
	assert(start_rank >= 0 &&
		   start_rank < parops->g_np);
	while (rc_seq_start_idx >
		   g_seq_offsets[start_rank]+l_seq_counts[start_rank])
		++start_rank;
	assert(start_rank >= 0 &&
		   start_rank < parops->g_np);

	uint64_t rc_seq_count_needed = avg_grid_seq_count;
	if (grid_rc_id == grid_rc_procs_count-1) // procs at last grid dim
		rc_seq_count_needed = g_seq_count-(avg_grid_seq_count*grid_rc_id);

	int			nbr_rank;
  	uint64_t	nbr_seq_start_idx, nbr_seq_end_idx;
  	uint64_t	count		  = 0;
  	uint64_t	seq_start_idx = rc_seq_start_idx;
	while (count < rc_seq_count_needed)
	{
		nbr_rank = start_rank;
		nbr_seq_start_idx = seq_start_idx - g_seq_offsets[start_rank];
		uint64_t remaining_needed = rc_seq_count_needed - count;
		uint64_t remaining_in_start_rank =
			l_seq_counts[start_rank] - nbr_seq_start_idx;

		if (remaining_needed >= remaining_in_start_rank)
		{
			count += remaining_in_start_rank;
			nbr_seq_end_idx = (seq_start_idx+remaining_in_start_rank-1) -
				g_seq_offsets[start_rank];
		}
		else
		{
			count += remaining_needed;
			nbr_seq_end_idx = (seq_start_idx+remaining_needed-1) -
				g_seq_offsets[start_rank];
		}

		my_nbrs.emplace_back(rc_flag, parops->g_rank, nbr_rank, 1001,
							 nbr_seq_start_idx, nbr_seq_end_idx);
		++start_rank;
		if (start_rank >= parops->g_np)
		{
			assert(count == rc_seq_count_needed);
			break;
		}

		seq_start_idx = g_seq_offsets[start_rank];
	}

	assert(count == rc_seq_count_needed);
}



void
DistFastaData::wait ()
{
	parops->tp->start_timer("dfd:MPI_Waitall(seqs)");
	MPI_Waitall(recv_nbrs_count, recv_nbrs_buffs_reqs, recv_nbrs_buffs_stats);
	MPI_Waitall(to_nbrs_count, to_nbrs_buffs_reqs, to_nbrs_buffs_stat);
	parops->tp->stop_timer("dfd:MPI_Waitall(seqs)");

	parops->tp->stop_timer("dfd:extract_recv_seqs");

	int recv_nbr_idx = 0;
	uint64_t rseq_cnt = 0, cseq_cnt = 0;
	for (auto &nbr : my_nbrs)
	{
		uint64_t nbr_seqs_count =
			(nbr.nbr_seq_end_idx - nbr.nbr_seq_start_idx) + 1;		
		if (nbr.rc_flag == 1)
			rseq_cnt += nbr_seqs_count;
		else
			cseq_cnt += nbr_seqs_count;
		
		if (nbr.nbr_rank == parops->g_rank) // already have meself
			continue;

		uint64_t recv_nbr_l_end = recv_nbrs_buff_lengths[recv_nbr_idx]-1;
		FastaData *recv_fd = new FastaData(recv_nbrs_buffs[recv_nbr_idx], k, 0,
										   recv_nbr_l_end);
		recv_fds_.push_back(recv_fd);
		++recv_nbr_idx;			
	}

	parops->tp->stop_timer("dfd:extract_recv_seqs");

	parops->logger->log("computed #row seqs " + to_string(rseq_cnt) +
						" #col seqs " + to_string(cseq_cnt));
	assert(rseq_cnt == (row_seq_end_idx-row_seq_start_idx)+1 &&
		   cseq_cnt == (col_seq_end_idx-col_seq_start_idx)+1);

	ready = true;

	if (parops->g_rank == 0)
		std::cout << "all seqs received, wait complete." << std::endl;
}



void
DistFastaData::collect_seqs
(
 	int br,
 	int bc
)
{
	string s_tmp;
	vector<pair<uint64_t, uint64_t>> roffsets =
		cbmat_offsets(g_seq_count, br, false);
	vector<pair<uint64_t, uint64_t>> coffsets =
		cbmat_offsets(g_seq_count, bc, true);

	bl_rseq_offsets_.resize(br+1);
	bl_rseq_offsets_[0] = 0;
	for (int i = 0; i < br; ++i)
		bl_rseq_offsets_[i+1] = bl_rseq_offsets_[i]+roffsets[i].second;

	bl_cseq_offsets_.resize(bc+1);
	bl_cseq_offsets_[0] = 0;
	for (int i = 0; i < bc; ++i)
		bl_cseq_offsets_[i+1] = bl_cseq_offsets_[i]+coffsets[i].second;

	s_tmp = "";
	for (auto el : bl_rseq_offsets_)
		s_tmp += to_string(el) + " ";
	parops->logger->log("block row offsets: " + s_tmp);

	s_tmp = "";
	for (auto el : bl_cseq_offsets_)
		s_tmp += to_string(el) + " ";
	parops->logger->log("block col offsets: " + s_tmp);

	
	// form comm neighbors
	int tag_idx = 0;
	s_tmp = "";
	for (int bi = 0; bi<br+bc; ++bi)
	{
		bool	 is_row_block = (bi<br) ? true : false;
		int		 cur_block	  = (bi<br) ? bi : bi-br;
		auto	&b_offsets	  = is_row_block ? roffsets : coffsets;

		// cannot avoid  diag proc since block dims can differ
		// if (!is_row_block &&
		// 	(parops->grid->GetRankInProcCol() ==
		// 	 parops->grid->GetRankInProcRow())) // diag proc
		// 	continue;

		uint64_t seq_idx_beg = b_offsets[cur_block].first;
		uint64_t seq_idx_end = seq_idx_beg + b_offsets[cur_block].second; // exc

		s_tmp += "\n" + (is_row_block ? string("row block " ) :
						 string("col block ")) +
			to_string(cur_block) + " seq range [" + to_string(seq_idx_beg) +
			", " + to_string(seq_idx_end) + ")";
		
		// search for the range
		assert(seq_idx_beg >= 0 && seq_idx_beg < g_seq_count);
		auto iter = upper_bound(g_seq_offsets, g_seq_offsets+parops->g_np+1,
								seq_idx_beg);
		int beg_proc = iter-g_seq_offsets-1;
		s_tmp += " beg proc " + to_string(beg_proc);

		assert(seq_idx_end < g_seq_count-1);
		iter = upper_bound(g_seq_offsets, g_seq_offsets+parops->g_np+1,
						   seq_idx_end-1); // inc
		int end_proc = iter-g_seq_offsets-1; // inc
		s_tmp += " end proc " + to_string(end_proc);

		// add comm neighbors
		uint64_t cur_beg = seq_idx_beg;
		s_tmp += " ||| ";
		for (int cur_proc = beg_proc; cur_proc<=end_proc-1; ++cur_proc)
		{
			bl_nbrs_.emplace_back(is_row_block, parops->g_rank, cur_proc,
								  tag_idx++,
								  cur_beg-g_seq_offsets[cur_proc],
								  g_seq_offsets[cur_proc+1]-
								  g_seq_offsets[cur_proc]-1); // inc
			s_tmp += "[" + to_string(cur_beg) + ", " +
				to_string(g_seq_offsets[cur_proc+1]-1) + "] from " +
				to_string(cur_proc) + " - ";
			cur_beg = g_seq_offsets[cur_proc+1];			
		}

		bl_nbrs_.emplace_back(is_row_block, parops->g_rank, end_proc, tag_idx++,
							  cur_beg-g_seq_offsets[end_proc],
							  seq_idx_end-g_seq_offsets[end_proc]-1); // inc
		s_tmp += "[" + to_string(cur_beg) + ", " +
				to_string(seq_idx_end-1) + "] from " +
				to_string(end_proc) + " - ";		
	}

	parops->logger->log(s_tmp);		


	////////////////////////////////////////////////////////////////////////////

	// communicate neighbor and size data
	int block_length[6] = {1, 1, 1, 1, 1, 1};
	MPI_Aint displacement[6] = {offsetof(NbrData, rc_flag),
								offsetof(NbrData, owner_rank),
								offsetof(NbrData, nbr_rank),
								offsetof(NbrData, tag),
								offsetof(NbrData, nbr_seq_start_idx),
								offsetof(NbrData, nbr_seq_end_idx)};
	MPI_Datatype types[] = {MPI_UNSIGNED_SHORT, MPI_INT, MPI_INT, MPI_INT,
							MPI_UINT64_T, MPI_UINT64_T};
	MPI_Datatype MPI_NbrData;
	MPI_Type_create_struct(6, block_length, displacement, types, &MPI_NbrData);
	MPI_Type_commit(&MPI_NbrData);

	for (auto &nbr : bl_nbrs_)
	{
		if (nbr.nbr_rank != parops->g_rank)
			bl_cnbrs_.emplace_back(&nbr, 0);
	}

	// bl_nrseqs_ = 0;
	// bl_ncseqs_ = 0;
	// for_each(bl_nbrs_.begin(), bl_nbrs_.end(),
	// 		 [&]
	// 		 (const NbrData &nbr)
	// 		 {
	// 			 uint64_t tmp =
	// 				 nbr.nbr_seq_end_idx - nbr.nbr_seq_start_idx + 1;
	// 			 if (nbr.rc_flag == 1)
	// 				 bl_nrseqs_ += tmp;
	// 			 else
	// 				 bl_ncseqs_ += tmp;
	// 		 });

	parops->logger->log("number of row seqs " +
						to_string(bl_rseq_offsets_[bl_rseq_offsets_.size()-1]) +
						" number of col seqs " +
						to_string(bl_cseq_offsets_[bl_cseq_offsets_.size()-1]));
	parops->logger->log("number of messages to be received " +
						to_string(bl_cnbrs_.size()));

	// sizes of the buffers for the seqs to be received
	MPI_Request *rbuflen_reqs = new MPI_Request[bl_cnbrs_.size()];
	for (int i = 0; i < bl_cnbrs_.size(); ++i)
		MPI_Irecv(&(bl_cnbrs_[i].second),
				  1,
				  MPI_UINT64_T,
				  bl_cnbrs_[i].first->nbr_rank,
				  1001+bl_cnbrs_[i].first->tag,
				  MPI_COMM_WORLD,
				  rbuflen_reqs + i);

	parops->logger->log("communicating neighbor metadata");

	// gather all neighbor metadata first
	int	 my_nbrs_count	 = bl_nbrs_.size();
	int *all_nbrs_counts = new int[parops->g_np];
	MPI_Allgather(&my_nbrs_count,
				  1,
				  MPI_INT,
				  all_nbrs_counts,
				  1,
				  MPI_INT,
				  MPI_COMM_WORLD);

	int		 all_nbrs_count	  = 0;
	auto	*all_nbrs_displas = new int[parops->g_np];
	all_nbrs_displas[0]		  = 0;
	for (int i = 0; i < parops->g_np; ++i)
	{
		all_nbrs_count += all_nbrs_counts[i];
		if (i > 0)
			all_nbrs_displas[i] = all_nbrs_displas[i-1] + all_nbrs_counts[i-1];
	}

	auto *all_nbrs = new NbrData[all_nbrs_count];
	MPI_Allgatherv(&bl_nbrs_[0],
				   my_nbrs_count,
				   MPI_NbrData,
				   all_nbrs,
				   all_nbrs_counts,
				   all_nbrs_displas,
				   MPI_NbrData,
				   MPI_COMM_WORLD);

	sort(all_nbrs, all_nbrs + all_nbrs_count,
		 [](const NbrData &a, const NbrData &b) -> bool
		 {
			 return a.nbr_rank < b.nbr_rank;
		 });

	vector<uint64_t> send_lengths, send_start_offsets;
	vector<NbrData *> to_nbrs;
	bl_nsnbrs_ = 0;
	s_tmp = "\n";
	for (int i = 0; i < all_nbrs_count; ++i)
	{
		NbrData &nbr = all_nbrs[i];
		if ((nbr.nbr_rank != parops->g_rank) || // doesn't request from me
			(nbr.owner_rank == parops->g_rank))	// local seq
			continue;

		uint64_t len, start_offset, end_offset;
		fd->buffer_size(nbr.nbr_seq_start_idx, nbr.nbr_seq_end_idx,
						len, start_offset, end_offset);
		send_lengths.push_back(len);
		send_start_offsets.push_back(start_offset);
		to_nbrs.push_back(&nbr);
		++bl_nsnbrs_;

		s_tmp += "to P " + to_string(nbr.owner_rank) +
			" [" + to_string(nbr.nbr_seq_start_idx) + ", " +
			to_string(nbr.nbr_seq_end_idx) +
			"] row? " + to_string(nbr.rc_flag) + 
			" size " +
			to_string(len/static_cast<double>(1<<20)) + "\n";
	}

	parops->logger->log("Seq send sizes (MB) " + s_tmp);

	MPI_Request *sbuflen_reqs = new MPI_Request[bl_nsnbrs_];
	for (int i = 0; i < bl_nsnbrs_; ++i)
		MPI_Isend(&send_lengths[i],
				  1,
				  MPI_UINT64_T,
				  to_nbrs[i]->owner_rank,
				  1001+to_nbrs[i]->tag,
				  MPI_COMM_WORLD,
				  sbuflen_reqs+i);

	// wait for size comm
	MPI_Status *rstatus = new MPI_Status[bl_cnbrs_.size()];
	MPI_Status *sstatus = new MPI_Status[bl_nsnbrs_];
	MPI_Waitall(bl_cnbrs_.size(), rbuflen_reqs, rstatus);
	MPI_Waitall(bl_nsnbrs_, sbuflen_reqs, sstatus);

	s_tmp = "";
	for (int i = 0; i < bl_cnbrs_.size(); ++i)
		s_tmp += "from P " + to_string(bl_cnbrs_[i].first->nbr_rank) +
			" " + to_string(bl_cnbrs_[i].second/
						  static_cast<double>(1<<20)) + " | ";
	parops->logger->log("Seq recv sizes (MB) " + s_tmp);

	uint64_t tot_ssz = 0, tot_rsz = 0;
	for (int i = 0; i < bl_nsnbrs_; ++i)
	{
		tot_ssz += send_lengths[i];
		if (send_lengths[i] > std::numeric_limits<int>::max())
			parops->logger->log("Send seq size larger than INT_MAX: " +
								to_string(send_lengths[i]) + " bytes from P " +
								to_string(parops->g_rank) + " to P " +
								to_string(to_nbrs[i]->owner_rank),
								Logger::LogLevel::WARNING);
	}

	for (int i = 0; i < bl_cnbrs_.size(); ++i)
	{
		tot_rsz += bl_cnbrs_[i].second;
		if (bl_cnbrs_[i].second > std::numeric_limits<int>::max())
			parops->logger->log("Recv seq size larger than INT_MAX: " +
								to_string(bl_cnbrs_[i].second) +
								" bytes from P " +
								to_string(bl_cnbrs_[i].first->nbr_rank) +
								" to P " +
								to_string(parops->g_rank),
								Logger::LogLevel::WARNING);
	}

	parops->logger->log("total send size (MB) " +
						to_string(tot_ssz/static_cast<double>(1<<20)));
	parops->logger->log("total recv size (MB) " +
						to_string(tot_rsz/static_cast<double>(1<<20)));

	// @OGUZ-TODO do multiple send/recvs below if buff length > int max

	parops->logger->log("starting communication of seqs");

	#if PASTIS_DBG_LVL > 0
	for (int i = 0; i < bl_cnbrs_.size(); ++i)
		parops->bytes_alloc += bl_cnbrs_[i].second * sizeof(char);
	#endif

	for (int i = 0; i < bl_cnbrs_.size(); ++i)
	{
		char *tmp = new char[bl_cnbrs_[i].second];
		bl_rbufs_.emplace_back(tmp, nullptr);
	}

	bl_rreqs_ = new MPI_Request[bl_cnbrs_.size()];
	for (int i = 0; i < bl_cnbrs_.size(); ++i)
		MPI_Irecv(bl_rbufs_[i].first,
				  static_cast<int>(bl_cnbrs_[i].second),
				  MPI_CHAR,
				  bl_cnbrs_[i].first->nbr_rank,
				  100001+bl_cnbrs_[i].first->tag,
				  MPI_COMM_WORLD,
				  bl_rreqs_+i);

	bl_sreqs_ = new MPI_Request[bl_nsnbrs_];
	for (int i = 0; i < bl_nsnbrs_; ++i)
		MPI_Isend(fd->buffer() + send_start_offsets[i],
				  static_cast<int>(send_lengths[i]),
				  MPI_CHAR,
				  to_nbrs[i]->owner_rank,
				  100001+to_nbrs[i]->tag,
				  MPI_COMM_WORLD,
				  bl_sreqs_+i);
	

	delete[] rbuflen_reqs;
	delete[] all_nbrs_counts;
	delete[] all_nbrs_displas;
	delete[] all_nbrs;
	delete[] sbuflen_reqs;
	delete[] rstatus;
	delete[] sstatus;
}



void
DistFastaData::bl_wait ()
{
	MPI_Status *bl_rstts = new MPI_Status[bl_cnbrs_.size()];
	MPI_Status *bl_sstts = new MPI_Status[bl_nsnbrs_];
	
	parops->tp->start_timer("dfd:MPI_Waitall(seqs)");
	MPI_Waitall(bl_cnbrs_.size(), bl_rreqs_, bl_rstts);
	MPI_Waitall(bl_nsnbrs_, bl_sreqs_, bl_sstts);
	parops->tp->stop_timer("dfd:MPI_Waitall(seqs)");

	parops->tp->start_timer("dfd:extract_recv_seqs");

	uint64_t rseq_cnt = 0, cseq_cnt = 0;
	int i = 0;
	for (auto &nbr : bl_nbrs_)
	{
		uint64_t nbr_seqs_count =
			(nbr.nbr_seq_end_idx-nbr.nbr_seq_start_idx)+1;		
		if (nbr.rc_flag == 1)
			rseq_cnt += nbr_seqs_count;
		else
			cseq_cnt += nbr_seqs_count;
		
		if (nbr.nbr_rank == parops->g_rank) // already have meself
			continue;

		uint64_t l_end = bl_cnbrs_[i].second-1;
		FastaData *tmp = new FastaData(bl_rbufs_[i].first, k, 0, l_end);
		bl_rbufs_[i].second = tmp;		
		++i;

		parops->logger->log("bl_wait nbr_seqs_count " + to_string(nbr_seqs_count) +
							" fd " + to_string(tmp->local_count()));
		assert(nbr_seqs_count == tmp->local_count());
	}

	parops->tp->stop_timer("dfd:extract_recv_seqs");

	parops->logger->log("computed #row seqs " + to_string(rseq_cnt) +
						" #col seqs " + to_string(cseq_cnt));
	assert(rseq_cnt == bl_rseq_offsets_[bl_rseq_offsets_.size()-1] &&
		   cseq_cnt == bl_cseq_offsets_[bl_cseq_offsets_.size()-1]);

	ready = true;

	if (parops->g_rank == 0)
		std::cout << "all seqs received, wait complete." << std::endl;

	delete[] bl_rstts;
	delete[] bl_sstts;
}

	
}
