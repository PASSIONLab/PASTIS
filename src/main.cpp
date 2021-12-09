#include <algorithm>
#include <cstdint>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "CombBLAS/CombBLAS.h"

#include "../inc/align/PWAlign.hpp"
#include "../inc/align/SeqanFullAligner.hpp"
#include "../inc/Alphabet.hpp"
#include "../inc/cxxopts.hpp"
#include "../inc/Constants.hpp"
#include "../inc/DistFastaData.hpp"
#include "../inc/DistPWRunner.hpp"
#include "../inc/kmer/Kmer.hpp"
#include "../inc/kmer/KmerOps.hpp"
#include "../inc/ParallelOps.hpp"
#include "../inc/sim_search.hpp"
#include "../inc/sr.hpp"
#include "../inc/TraceUtils.hpp"
#include "../inc/Types.hpp"
#include "../inc/util.hpp"



using std::cout;	using std::endl;	using std::string;	using std::ofstream;
using std::to_string;	using std::vector;	using std::shared_ptr;
using std::max;

extern shared_ptr<pastis::ParallelOps> parops; // global var



int
parse_args
(
    int					  argc,
	char				**argv,
	pastis::params_t	 &params
)
{
	cxxopts::Options options("PASTIS",
							 "A fast, distributed protein aligner");

	options.add_options()
    (pastis::CMD_OPTION_INPUT, pastis::CMD_OPTION_DESCRIPTION_INPUT,
     cxxopts::value<string>())
    (pastis::CMD_OPTION_INPUT_SEQ_COUNT, pastis::CMD_OPTION_DESCRIPTION_INPUT_SEQ_COUNT,
     cxxopts::value<uint64_t>())
    (pastis::CMD_OPTION_INPUT_OVERLAP, pastis::CMD_OPTION_DESCRIPTION_INPUT_OVERLAP,
     cxxopts::value<uint64_t>())
    (pastis::CMD_OPTION_SEED_COUNT, pastis::CMD_OPTION_DESCRIPTION_SEED_COUNT,
     cxxopts::value<int>())
    (pastis::CMD_OPTION_GAP_OPEN, pastis::CMD_OPTION_DESCRIPTION_GAP_OPEN,
     cxxopts::value<int>())
    (pastis::CMD_OPTION_GAP_EXT, pastis::CMD_OPTION_DESCRIPTION_GAP_EXT,
     cxxopts::value<int>())
    (pastis::CMD_OPTION_KMER_LENGTH, pastis::CMD_OPTION_DESCRIPTION_KMER_LENGTH,
     cxxopts::value<int>())
    (pastis::CMD_OPTION_KMER_STRIDE, pastis::CMD_OPTION_DESCRIPTION_KMER_STRID,
     cxxopts::value<int>())
    (pastis::CMD_OPTION_OVERLAP_FILE, pastis::CMD_OPTION_DESCRIPTION_OVERLAP_FILE,
     cxxopts::value<string>())
    (pastis::CMD_OPTION_ALIGN_FILE, pastis::CMD_OPTION_DESCRIPTION_ALIGN_FILE,
     cxxopts::value<string>())
    (pastis::CMD_OPTION_NO_ALIGN, pastis::CMD_OPTION_DESCRIPTION_NO_ALIGN)
    (pastis::CMD_OPTION_FULL_ALIGN, pastis::CMD_OPTION_DESCRIPTION_FULL_ALIGN)
    (pastis::CMD_OPTION_XDROP_ALIGN, pastis::CMD_OPTION_DESCRIPTION_XDROP_ALIGN,
     cxxopts::value<int>())
    (pastis::CMD_OPTION_BANDED_ALIGN, pastis::CMD_OPTION_DESCRIPTION_BANDED_ALIGN,
     cxxopts::value<int>())
	(pastis::CMD_OPTION_GPUBSW_ALIGN, pastis::CMD_OPTION_DESCRIPTION_GPUBSW_ALIGN)
    (pastis::CMD_OPTION_IDX_MAP, pastis::CMD_OPTION_DESCRIPTION_IDX_MAP,
     cxxopts::value<string>())
    (pastis::CMD_OPTION_ALPH, pastis::CMD_OPTION_DESCRIPTION_ALPH,
     cxxopts::value<string>())
    (pastis::CMD_OPTION_SUBS, pastis::CMD_OPTION_DESCRIPTION_SUBS,
     cxxopts::value<int>())
    (pastis::CMD_OPTION_AF_FREQ, pastis::CMD_OPTION_DESCRIPTION_AF_FREQ,
     cxxopts::value<int>())
	(pastis::CMD_OPTION_CKTHR, pastis::CMD_OPTION_DESCRIPTION_CKTHR,
     cxxopts::value<int>())
	(pastis::CMD_OPTION_MOSTHR, pastis::CMD_OPTION_DESCRIPTION_MOSTHR,
     cxxopts::value<float>())
	(pastis::CMD_OPTION_ALN_BATCH_SZ, pastis::CMD_OPTION_DESCRIPTION_ALN_BATCH_SZ,
     cxxopts::value<uint64_t>())
	(pastis::CMD_OPTION_SIM_BR, pastis::CMD_OPTION_DESCRIPTION_SIM_BR,
     cxxopts::value<int>())
    (pastis::CMD_OPTION_SIM_BC, pastis::CMD_OPTION_DESCRIPTION_SIM_BC,
     cxxopts::value<int>());

	// defaults
	params.write_overlaps	   = false;
	params.add_substitue_kmers = false;
	params.subk_count		   = 0;
	params.pw_aln			   = pastis::params_t::PwAln::ALN_NONE;
	params.aln_seqan_banded_hw = 5;	
	params.seed_count		   = 2;
	params.ckthr			   = 0;
	params.mosthr			   = -1.0f;
	params.idx_map_file 	   = "";
	params.afreq			   = 1e6;
	params.aln_batch_sz		   = 1e7;
	params.aln_cov_thr		   = 0.7;
	params.aln_ani_thr		   = 30;
	params.br				   = -1; // default not-blocked
	params.bc				   = -1;

	bool	is_world_rank0 = parops->g_rank == 0;
	auto	result		   = options.parse(argc, argv);

	if (result.count(pastis::CMD_OPTION_INPUT))
		params.input_file = result[pastis::CMD_OPTION_INPUT].as<string>();
	else
	{
		if (is_world_rank0)
			cout << "ERROR: Input file not specified." << endl;
		return -1;
    }  

    if (result.count(pastis::CMD_OPTION_IDX_MAP))
	{
        params.idx_map_file = result[pastis::CMD_OPTION_IDX_MAP].as<string>();
	}

  	if (result.count(pastis::CMD_OPTION_INPUT_SEQ_COUNT))
	{
		params.seq_count =
			result[pastis::CMD_OPTION_INPUT_SEQ_COUNT].as<uint64_t>();
	}
	else
	{
    	if (is_world_rank0)
			cout << "ERROR: Input sequence count not specified" << endl;
    	return -1;
  	}

  	if (result.count(pastis::CMD_OPTION_INPUT_OVERLAP))
		params.input_overlap =
			result[pastis::CMD_OPTION_INPUT_OVERLAP].as<uint64_t>();
	else
		params.input_overlap = 10000;

  	if (result.count(pastis::CMD_OPTION_SEED_COUNT))
		params.seed_count = result[pastis::CMD_OPTION_SEED_COUNT].as<int>();
	else
		params.seed_count = 2;

  	if (result.count(pastis::CMD_OPTION_GAP_OPEN))
		params.gap_open = result[pastis::CMD_OPTION_GAP_OPEN].as<int>();
	else
		params.gap_open = -11;

  	if (result.count(pastis::CMD_OPTION_GAP_EXT))
		params.gap_ext = result[pastis::CMD_OPTION_GAP_EXT].as<int>();
	else
		params.gap_ext = -2;

  	if (result.count(pastis::CMD_OPTION_KMER_LENGTH))
		params.klength = result[pastis::CMD_OPTION_KMER_LENGTH].as<int>();
	else
		params.klength = 6;

  	if (result.count(pastis::CMD_OPTION_KMER_STRIDE))
		params.kstride = result[pastis::CMD_OPTION_KMER_STRIDE].as<int>();
	else
		params.kstride = 1;

  	if (result.count(pastis::CMD_OPTION_OVERLAP_FILE))
	{
    	params.write_overlaps = true;
    	params.overlap_file =
			result[pastis::CMD_OPTION_OVERLAP_FILE].as<string>();
	}

  	if (result.count(pastis::CMD_OPTION_ALIGN_FILE))
		params.align_file = result[pastis::CMD_OPTION_ALIGN_FILE].as<string>();

  	if (result.count(pastis::CMD_OPTION_NO_ALIGN))
		params.pw_aln = pastis::params_t::PwAln::ALN_NONE;

  	if (result.count(pastis::CMD_OPTION_FULL_ALIGN))
		params.pw_aln = pastis::params_t::PwAln::ALN_SEQAN_FULL;

  	if (result.count(pastis::CMD_OPTION_BANDED_ALIGN))
	{
    	params.pw_aln = pastis::params_t::PwAln::ALN_SEQAN_BANDED;
    	params.aln_seqan_banded_hw =
			result[pastis::CMD_OPTION_BANDED_ALIGN].as<int>();
	}

  	if (result.count(pastis::CMD_OPTION_XDROP_ALIGN))
	{
    	params.pw_aln = pastis::params_t::PwAln::ALN_SEQAN_XDROP;
    	params.aln_seqan_xdrop =
			result[pastis::CMD_OPTION_XDROP_ALIGN].as<int>();
	}
  
  	if (result.count(pastis::CMD_OPTION_GPUBSW_ALIGN))
		params.pw_aln = pastis::params_t::PwAln::ALN_ADEPT_GPUBSW;

  	if (result.count(pastis::CMD_OPTION_ALPH))
	{
		params.alph_str = result[pastis::CMD_OPTION_ALPH].as<string>();
		if (params.alph_str == string("pdefault"))
			params.alph = new pastis::DefaultProteinAlph();
		else if (params.alph_str == string("murphy10"))
			params.alph = new pastis::Murphy10ProteinAlph();
		else if (params.alph_str == string("dssp10"))
			params.alph = new pastis::DSSP10ProteinAlph();
		else if (params.alph_str == string("gbmr10"))
			params.alph = new pastis::GBMR10ProteinAlph();
		else if (params.alph_str == string("td10"))
			params.alph = new pastis::TD10ProteinAlph();
		else if (params.alph_str == string("diamond"))
			params.alph = new pastis::DiamondProteinAlph();
		else
		{
		  	if (is_world_rank0)
				cout << "ERROR: Unsupported alphabet type "
					 << params.alph_str << endl;
		}
	}
	else
	{
		params.alph = new pastis::DefaultProteinAlph();
		params.alph_str = "pdefault";
	}


  	if (result.count(pastis::CMD_OPTION_SUBS))
	{
    	params.add_substitue_kmers = true;
    	params.subk_count = result[pastis::CMD_OPTION_SUBS].as<int>();
	}

	if (result.count(pastis::CMD_OPTION_AF_FREQ))
		params.afreq = result[pastis::CMD_OPTION_AF_FREQ].as<int>();

	if (result.count(pastis::CMD_OPTION_CKTHR))
		params.ckthr = result[pastis::CMD_OPTION_CKTHR].as<int>();

	if (result.count(pastis::CMD_OPTION_MOSTHR))
		params.mosthr = result[pastis::CMD_OPTION_MOSTHR].as<float>();

	if (result.count(pastis::CMD_OPTION_ALN_BATCH_SZ))
		params.aln_batch_sz =
			result[pastis::CMD_OPTION_ALN_BATCH_SZ].as<uint64_t>();

	if (result.count(pastis::CMD_OPTION_SIM_BR))
		params.br = result[pastis::CMD_OPTION_SIM_BR].as<int>();

	if (result.count(pastis::CMD_OPTION_SIM_BC))
		params.bc = result[pastis::CMD_OPTION_SIM_BC].as<int>();	

  	return 0;
}



void
params_to_str
(
 	const pastis::params_t	&params,
    string					&s
	
)
{
	vector<string> params_str = {
		"Input file (-i)",
		"Original sequence count (-c)",
		"Kmer length (k)",
		"Kmer stride (s)",
		"Overlap in bytes (-O)",
		"Max seed count (--sc)",
		"Gap open penalty (-g)",
		"Gap extension penalty (-e)",
		"Overlap file (--of)",
		"Alignment file (--af)",
		"Alignment write frequency (--afreq)",
		"Pairwise seq aligner (--na | --sfa | --sxa | --sba | -absw)",
		"Index map (--idxmap)",
		"Alphabet (--alph)",
		"Use substitute kmers (--subs)",
		"Common k-mer threshold (--ckthr)",
		"Max overlap score threshold (--mosthr)",
		"Batch alignment size (--bsz)",
		"Block mult row dim (--br)",
		"Block mult col dim (--bc)"
	};

	auto bool_to_str = [] (bool b)
		{
			return b ? "True" : "False";
		};

	string aln_str = "";
	if (params.pw_aln == pastis::params_t::PwAln::ALN_NONE)
		aln_str = "No alignment";
	else if (params.pw_aln == pastis::params_t::PwAln::ALN_SEQAN_FULL)
		aln_str = "Seqan full";
	else if (params.pw_aln == pastis::params_t::PwAln::ALN_SEQAN_XDROP)
		aln_str = "Seqan xdrop (xdrop value " +
			to_string(params.aln_seqan_xdrop) + ")";
	else if (params.pw_aln == pastis::params_t::PwAln::ALN_SEQAN_BANDED)
		aln_str = "Seqan banded (half band width " +
			to_string(params.aln_seqan_banded_hw) + ")";
	else if (params.pw_aln == pastis::params_t::PwAln::ALN_ADEPT_GPUBSW)
		aln_str = "ADEPT BSW (GPU)";

	vector<string> vals = {
		params.input_file,
		to_string(params.seq_count),
		to_string(params.klength),
		to_string(params.kstride),
		to_string(params.input_overlap),
		to_string(params.seed_count),
		to_string(params.gap_open),
		to_string(params.gap_ext),
		!params.overlap_file.empty() ? params.overlap_file : "None",
		!params.align_file.empty() ? params.align_file : "None",
		!params.align_file.empty() ? to_string(params.afreq) : "None",
		aln_str,
		!params.idx_map_file.empty() ? params.idx_map_file : "None",
		params.alph_str,
		bool_to_str(params.add_substitue_kmers) + (params.add_substitue_kmers ?
											" | sub kmers: " +
											to_string(params.subk_count) : ""),
		to_string(params.ckthr),
		to_string(params.mosthr),
		to_string(params.aln_batch_sz),
		to_string(params.br),
		to_string(params.bc)
	  };

	int max_length = 0;
	for (const auto &param : params_str)
		max_length = max(static_cast<int>(param.size()), max_length);

	s.append("Parameters...\n");
	for (size_t i = 0; i < params_str.size(); ++i)
	{
		string param = params_str[i];
		s.append("  ").append(param).append(": ")
			.append(string(max_length-param.size(), ' '))
			.append(vals[i]).append("\n");
	}  
}
	


int
main
(
    int		  argc,
	char	**argv
)
{
	parops = pastis::ParallelOps::init(&argc, &argv);
	bool is_log_active;
	string s_tmp;
	
	// process arguments	
	pastis::params_t params;
	int ret = parse_args(argc, argv, params);
	if (ret < 0)
	{
	  	parops->teardown_parallelism();
	  	return ret;
	}

	int nthreads = 1;
  	#ifdef THREADED
  	#pragma omp parallel
	{
		nthreads = omp_get_num_threads();
	}
  	#endif

	if (parops->g_rank == 0)
		cout << "Process Grid (p x p x t): "
			 << sqrt(parops->g_np) << " x " << sqrt(parops->g_np) << " x "
			 << nthreads << endl;

	// logger and timer
	is_log_active = true;
	// if (parops->g_rank == 0)
	// 	is_log_active = true;
	s_tmp = to_string(parops->g_rank);
	parops->logger = pastis::Logger::instantiate
		(is_log_active,
		 "pastis-" + string(4-s_tmp.size(), '0') + s_tmp + ".log");

	parops->tp->start_timer("main");
	s_tmp = "\nINFO: Program started on ";
	s_tmp.append(parops->tp->get_time_str("start_main")).append("\n");
	params_to_str(params, s_tmp);
	if (parops->g_rank == 0)
		cout << s_tmp << endl;

	if (params.pw_aln == pastis::params_t::PwAln::ALN_NONE ||
		params.pw_aln == pastis::params_t::PwAln::ALN_SEQAN_FULL ||
		params.pw_aln == pastis::params_t::PwAln::ALN_ADEPT_GPUBSW)
		pastis::compute_sim_mat<pastis::CommonKmerLight>(params);
	else if (params.pw_aln == pastis::params_t::PwAln::ALN_SEQAN_XDROP)
		pastis::compute_sim_mat<pastis::CommonKmerLoc>(params);	

	parops->tp->stop_timer("main");
	if (parops->g_rank == 0)
		cout << parops->tp->to_string() << endl;
	parops->logger->log("\n" + parops->tp->to_string());
	parops->teardown_parallelism();
	


	return (EXIT_SUCCESS);
}


