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
  	string names[26] =
		{
			"total",
			"fasta",
			"fasta|io-seqs",
			"fasta|local-setup",
			"fasta|seq-comm-setup",
			"seq-kmer-mat",
			"seq-kmer-mat|add-kmers",
			"seq-kmer-mat|assemble",
			"sparse",
			"tr",
			"sim-search",
			"seq-comm-wait",
			"block-split",			
			"construct-seqs",
			"mult-align",
			"mat-mult",
			"prune",
			"mat-formation",
			"sim-mat-assemble",			
			"sim-mat-io",			
			"align",
			"align|form-tuples",
			"align|batch",
			"align|pre",
			"align|cpu",
			"align|post"
		};


	// block nnz stats added cumulatively - entire matrix info
	string mat_stat_names[6] =
		{
		 "c_nnz",
		 "prune_sym",
		 "kmer_thr",
		 "sim_thr",
		 "aln_pairs",			// per-process
		 "aln_pair_lens"		// per-process
		};
	unordered_map<string, double> mat_stats;


	
	TimePod ()
	{
		for (string &s : mat_stat_names)
			mat_stats[s] = 0.0;
	}



	void
	get_proc_stat (double val, double &val_min, double &val_max,
				   double &val_avg, double &limb)
	{
		int np;
		MPI_Comm_size(MPI_COMM_WORLD, &np);
		MPI_Reduce(&val, &val_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
		MPI_Reduce(&val, &val_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		double tmp;
		MPI_Reduce(&val, &tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		val_avg = tmp / np;
		limb = 0;
		if (val_avg != 0)
			limb = 100 * ((val_max - val_avg) / val_avg);
	}
	

	
  	string
  	to_string ()
	{
		double val_min, val_max, val_avg, val_imb;

		// string str = "\nINFO: Program timings ...\n";
		string str(80, '*');
		str.append("\n");
		str.append("Timing in msec\n");
		for (const auto &name : names)
		{			
			int tmp = std::count(name.begin(), name.end(), '|');
			str.append((tmp+1)*2, ' ');
			double tmp2 = 0.0;
			if (elapsed.find(name) != elapsed.end())
				tmp2 = elapsed[name].second;
			
			get_proc_stat(tmp2, val_min, val_max, val_avg, val_imb);
			std::stringstream ss;
			ss << " --- "
			   << std::fixed << std::setprecision(3)
			   << val_min 
			   << " "
			   << val_max 
			   << " "
			   << val_avg 
			   << " "
			   << std::fixed << std::setprecision(2)
			   << val_imb ;
			str.append(name).append(" ")
				.append(std::to_string(tmp2))
				.append(ss.str());
			
			str.append("\n");
		}

		for (auto &s : mat_stat_names)
		{
			std::stringstream ss;
			ss << "  " << s << " "
			   << std::fixed << std::setprecision(0)
			   << mat_stats[s];
			if (s == "aln_pairs" || s == "aln_pair_lens")
			{
				get_proc_stat(mat_stats[s], val_min, val_max, val_avg, val_imb);
				
				ss << " --- "
				   << std::fixed << std::setprecision(0)
				   << val_min 
				   << " "
				   << val_max 
				   << " "
				   << val_avg 
				   << " "
				   << std::fixed << std::setprecision(2)
				   << val_imb;				
			}
			str.append(ss.str());
			str.append("\n");			
		}

		str += string(80, '*');
		str.append("\n");
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



	void
	add_to_timer (const string &s, double arg)
	{
		auto tmp = std::chrono::system_clock::now();
		if (elapsed.find(s) == elapsed.end())
			elapsed[s] = pair<ticks_t, double>(tmp, 0.0);
		elapsed[s].second += arg;
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
