// Created by Saliya Ekanayake on 2/10/19.

#pragma once


#include <chrono>
#include <ratio>
#include <string>
#include <unordered_map>
#include <utility>

using std::string;	using std::unordered_map;	using std::pair;

typedef std::chrono::duration<double, std::milli> ms_t;
typedef std::chrono::duration<double, std::micro> micros_t;
typedef std::chrono::time_point<std::chrono::high_resolution_clock> ticks_t;
typedef std::chrono::high_resolution_clock hrc_t;



namespace
pastis
{

struct
TimePod
{
  	// unordered_map<string, ticks_t> times;
	// unordered_map<string, double> elapsed;
	unordered_map<string, pair<ticks_t, double>> elapsed;
  	string names[29] =
		{
			"main",
			"main:newDFD",
			"dfd:pfr->read_fasta",
			"dfd:new_FD",
			"kmerop:gen_A:loop_add_kmers",
			"kmerop:gen_A:spMatA",
			"main:genA",
			"main:AT",
			"kmerop:gen_S:find_sub_kmers",
			"kmerop:gen_S:spMatS",
			"main:genS",
			"main:AxS",
			"main:(AS)AT",
			"main:dfd->wait",
			"dfd:MPI_Waitall(seqs)",
			"dfd:extract_recv_seqs",
			"main:sim->align",
			"main:sim->mult_align",
			"main:sim->write_overlaps",
			"sim:construct_seqs",			
			"sim:prune",
			"sim:align_all",
			"sim:align_pre",
			"sim:align",
			"sim:align_post",
			"sim:mat_formation",
			"sim:mat_mul",
			"sim:block_split",
			"sim:block_all"
		};


	
  	string
  	to_string ()
	{
		// string str = "\nINFO: Program timings ...\n";
		// ticks_t ts, te;
		// for (const auto &name : names)
		// {
		// 	if (times.find("start_" + name) != times.end())
		// 	{
		// 	  	ts = times["start_" + name];
		// 	  	te = times["end_" + name];
		// 	  	str.append("  ").append(name).append(":")
		// 			.append(std::to_string((ms_t(te - ts)).count()))
		// 			.append(" ms\n");
		// 	}
		// 	else
		// 		str.append("  ").append(name).append(" not evaluated.\n");
		// }
		// return str;

		string str = "\nINFO: Program timings ...\n";
		for (const auto &name : names)
		{
			if (elapsed.find(name) != elapsed.end())
				str.append("  ").append(name).append(":")
					.append(std::to_string(elapsed[name].second))
					.append(" ms\n");
			else
				str.append("  ").append(name).append(" not evaluated.\n");
		}
		return str;
	}



	void
	start_timer (const string &s)
	{
		auto tmp = std::chrono::system_clock::now();
		if (elapsed.find(s) == elapsed.end())
			elapsed[s] = pair<ticks_t, double>(tmp, 0.0);
		else
			elapsed[s].first = tmp;
	}



	void
	stop_timer (const string &s)
	{
		auto tmp = std::chrono::system_clock::now();
		if (elapsed.find(s) != elapsed.end())
			elapsed[s].second += ms_t(tmp-elapsed[s].first).count();
	}


	
	// void
	// add_time_point (const string &s)
	// {
	// 	times[s] = std::chrono::system_clock::now();
	// }



	string
	get_time_str (const string &s)
	{
		std::time_t tmp = std::chrono::system_clock::to_time_t(elapsed[s].first);
		return std::ctime(&tmp);
	}
};

}
