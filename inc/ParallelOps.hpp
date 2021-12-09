//
// Created by Saliya Ekanayake on 12/17/18.
//

#pragma once

#include <chrono>
#include <memory>
#include <string>

#include "mpi.h"
#include "CombBLAS/CombBLAS.h"
#include "CombBLAS/CommGrid.h"

#include "../inc/logger.hpp"
#include "../inc/TraceUtils.hpp"


using std::shared_ptr;


namespace
pastis
{

class
ParallelOps
{

private:
	
	static shared_ptr<ParallelOps> instance;


	ParallelOps (int world_proc_rank, int world_procs_count);




public:

	int								g_rank;
  	int								g_np;
  	shared_ptr<combblas::CommGrid>	grid;
	shared_ptr<Logger>				logger;
	uint64_t						bytes_alloc;
	shared_ptr<TimePod>				tp;



	~ParallelOps ();



	static const
	shared_ptr<ParallelOps>
	init (int *argc, char ***argv);

	void
	teardown_parallelism ();

	void
	write_file_in_parallel (const char *file, const std::string &local_data);



	void
	info (const std::string &s)
	{
		if (g_rank == 0)
			std::cout << s << std::endl;
	}
	
};

}
