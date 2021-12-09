// Created by Saliya Ekanayake on 12/17/18.

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

#include "mpi.h"

#include "../inc/ParallelOps.hpp"

using std::accumulate;	using std::string;



namespace
pastis
{

std::shared_ptr<ParallelOps> ParallelOps::instance = nullptr;



const
std::shared_ptr<ParallelOps>
ParallelOps::init
(
    int		  *argc,
	char	***argv
)
{
	if (instance != nullptr)
		return instance;

	int rank, count;
	MPI_Init(argc, argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &count);
	instance = std::shared_ptr<ParallelOps>(new ParallelOps(rank, count));
	return instance;
}



void
ParallelOps::teardown_parallelism ()
{
  	int flag;
  	MPI_Initialized(&flag);
  	if (!flag)
		MPI_Finalize();
}



ParallelOps::ParallelOps
(
    int world_proc_rank,
	int world_procs_count
) :
	g_rank(world_proc_rank),
	g_np(world_procs_count)
{

  auto	grid_size = static_cast<int>(sqrt(world_procs_count));
  grid = std::make_shared<combblas::CommGrid>(MPI_COMM_WORLD,
											  grid_size, grid_size);
  bytes_alloc = 0;
  tp = std::make_shared<TimePod>();
}


	
ParallelOps::~ParallelOps ()
{
    int flag;
	MPI_Initialized(&flag);
	if (!flag)
		MPI_Finalize();
}



void
ParallelOps::write_file_in_parallel
(
     const char     *file,
	 const string	&local_data
)
{
	/*! The following code is adopted from CombBLAS at
     * https://people.eecs.berkeley.edu/~aydin/CombBLAS/html/_fully_dist_sp_vec_8cpp_source.html#l01310.
     */

	auto *bytes = new int64_t[g_np];
	bytes[g_rank] = local_data.size();
	MPI_Allgather(MPI_IN_PLACE, 1, combblas::MPIType<int64_t>(), bytes, 1,
                  combblas::MPIType<int64_t>(), MPI_COMM_WORLD);

	int64_t bytesuntil = accumulate(bytes, bytes+g_rank,
									static_cast<int64_t>(0));
	int64_t bytestotal = accumulate(bytes, bytes+g_np, static_cast<int64_t>(0));

	if (g_rank == 0)
	{
        // only leader rights the original file with no content
        std::ofstream ofs(file, std::ios::out);
        ofs.seekp(bytestotal - 1);
        // this will likely create a sparse file so the actual disks won't spin yet
        ofs.write("", 1);
        ofs.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);

    struct stat st;     // get file row_size
    if (stat(file, &st) == -1)
	{
        MPI_Abort(MPI_COMM_WORLD, NOFILE);
    }
    
    // Then everyone fills it
    FILE *ffinal;
    if ((ffinal = fopen(file, "r+")) == NULL)
	{
        printf("ERROR: Output file %s failed to open at process %d\n",
               file, g_rank);
        MPI_Abort(MPI_COMM_WORLD, NOFILE);
    }

    fseek(ffinal, bytesuntil, SEEK_SET);
    fwrite(local_data.c_str(), 1, bytes[g_rank], ffinal);
    fflush(ffinal);
    fclose(ffinal);
    delete[] bytes;	
}
	
}


std::shared_ptr<pastis::ParallelOps> parops = nullptr; // global var
