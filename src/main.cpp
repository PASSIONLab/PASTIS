#include <iostream>
#include <cmath>
// #include <boost/uuid/uuid.hpp>
// #include <boost/uuid/uuid_generators.hpp>
// #include <boost/uuid/uuid_io.hpp>
#include "../include/Constants.hpp"
#include "../include/ParallelOps.hpp"
#include "../include/ParallelFastaReader.hpp"
#include "../include/Alphabet.hpp"
#include "../include/Utils.hpp"
#include "../include/DistributedPairwiseRunner.hpp"
#include "CombBLAS/CombBLAS.h"
#include "../include/cxxopts.hpp"
#include "../include/pw/SeedExtendXdrop.hpp"
#include "seqan/score/score_matrix_data.h"
#include "../include/pw/OverlapFinder.hpp"
#include "../include/pw/FullAligner.hpp"
#include "../include/pw/BandedAligner.hpp"
#include "../include/kmer/KmerOps.hpp"
#include "../include/kmer/KmerIntersectSR.hpp"
#include "../include/kmer/SubKmerIntersectSR.hpp"
#include "../include/kmer/sr.hpp"
#include <map>
#include <fstream>

/*! Namespace declarations */
using namespace combblas;

/*! Type definitions */
//typedef KmerIntersect<ushort, CommonKmers> KmerIntersectSR_t;
typedef pastis::KmerIntersect<pastis::MatrixEntry, pastis::CommonKmers> KmerIntersectSR_t;
typedef pastis::SubKmerIntersect<pastis::MatrixEntry, pastis::MatrixEntry> SubKmerIntersectSR_t;

/*! Function signatures */
int parse_args(int argc, char **argv);

void pretty_print_config(std::string &append_to);

std::string get_padding(ushort count, std::string prefix);

/*! Global variables */
std::shared_ptr<ParallelOps> parops;
std::string input_file;
uint64_t input_overlap;
uint64_t seq_count;
int xdrop;
int gap_open;
int gap_ext;
ushort klength;
ushort kstride;

/*! Parameters related to outputting k-mer overlaps */
bool write_overlaps = false;
std::string overlap_file;
bool add_substitue_kmers = false;
int subk_count = 0;

/*! Parameters related to outputting alignment info */
std::string align_file;
int afreq;

/*! Don't perform alignments if this flag is set */
bool no_align = false;

/*! Perform full alignment */
bool full_align = false;

/*! Perform xdrop alignment */
bool xdrop_align = false;

/*! Perform banded alignment */
bool banded_align = false;
int banded_half_width = 5;

/*! File path to output global sequence index to original global sequence
 * index mapping */
std::string idx_map_file;

/*! Alphabet to use. */
Alphabet::type alph_t;

bool is_print_rank = false;
std::string print_str;

/*! Maximum number of common k-mers to keep */
int seed_count = 2;

/*! Logging information */
std::string job_name = "pastis";
std::string proc_log_file;
std::ofstream proc_log_stream;
int log_freq;

// Common k-mer threshold
int ckthr = 0;

// Score threshold
float mosthr = -1.0;


int
main
(
    int argc,
	char **argv
)
{
	parops = ParallelOps::init(&argc, &argv);
	int ret = parse_args(argc, argv);
	if (ret < 0)
	{
	  	parops->teardown_parallelism();
	  	return ret;
	}

	/*! bcast job id */
	// Sender
	int job_name_length = job_name.length();

	if (parops->world_proc_rank == 0)
	{
	  	MPI_Bcast(&job_name[0], job_name_length, MPI_CHAR, 0, MPI_COMM_WORLD);
	} else
	{
		char* buf = new char[job_name_length];
		MPI_Bcast(buf, job_name_length, MPI_CHAR, 0, MPI_COMM_WORLD);
		std::string s(buf);
		job_name = s;
		delete [] buf;
	}

	int nthreads = 1;
  	#ifdef THREADED
  	#pragma omp parallel
	{
		nthreads = omp_get_num_threads();
	}
  	#endif

	int nprocs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	if (myrank == 0)
	{
		std::cout << "Process Grid (p x p x t): "
				  << sqrt(nprocs) << " x " << sqrt(nprocs) << " x "
				  << nthreads << std::endl;
	}

	proc_log_file = job_name + "_rank_" +
		std::to_string(parops->world_proc_rank) + "_log.txt";
	proc_log_stream.open(proc_log_file);

	is_print_rank = (parops->world_proc_rank == 0);
	std::shared_ptr<TimePod> tp = std::make_shared<TimePod>();
	TraceUtils tu(is_print_rank);

	/* Print start time information */
	tp->times["start_main"] = std::chrono::system_clock::now();
	std::time_t start_prog_time =
		std::chrono::system_clock::to_time_t(tp->times["start_main"]);
	print_str = "\nINFO: Program started on ";
	print_str.append(std::ctime(&start_prog_time));
	print_str.append("\nINFO: Job ID ").append(job_name).append("\n");
	pretty_print_config(print_str);
	tu.print_str(print_str);

	/* Read and distribute fasta data */
	tp->times["start_main:newDFD()"] = std::chrono::system_clock::now();
	std::shared_ptr<DistributedFastaData> dfd
	  = std::make_shared<DistributedFastaData>(
		input_file.c_str(), idx_map_file.c_str(), input_overlap,
		klength, parops, tp, tu);
	tp->times["end_main:newDFD()"] = std::chrono::system_clock::now();

  	#ifndef NDEBUG
		//  TraceUtils::print_fasta_data(fd, parops);
  	#endif

	if (dfd->global_count() != seq_count)
	{
		uint64_t final_seq_count = dfd->global_count();
		print_str = "\nINFO: Modfied sequence count\n";
		print_str.append("  Final sequence count: ")
		  .append(std::to_string(final_seq_count))
		  .append(" (").append(
			std::to_string((((seq_count-final_seq_count) * 100.0) / seq_count)))
		  .append("% removed)");

		seq_count = dfd->global_count();
		print_str += "\n";
		tu.print_str(print_str);
	}

	/*! Create alphabet */
	Alphabet alph(alph_t);

	/*! Generate sequences by kmers matrix */
	std::unordered_set<pastis::Kmer, pastis::Kmer> local_kmers;
	tp->times["start_main:genA()"] = std::chrono::system_clock::now();
	PSpMat<pastis::MatrixEntry>::MPI_DCCols A =
		pastis::KmerOps::generate_A(
			seq_count,dfd, klength, kstride,
			alph, parops, tp, local_kmers);
	tu.print_str("Matrix A: ");
	tu.print_str("\nLoad imbalance: " +
				 std::to_string(A.LoadImbalance()) + "\n");
	tp->times["end_main:genA()"] = std::chrono::system_clock::now();
	A.PrintInfo();

	// Transpose
	auto At = A;
	tp->times["start_main:At()"] = tp->times["end_main:genA()"];
	At.Transpose();
	tu.print_str("Matrix At: ");
	At.PrintInfo();
	tp->times["end_main:At()"] = std::chrono::system_clock::now();

	// Substitute matrix
	PSpMat<pastis::MatrixEntry>::MPI_DCCols* S = nullptr;
	if (add_substitue_kmers)
	{
	  	tp->times["start_main:genS()"] = std::chrono::system_clock::now();
	  	S = new PSpMat<pastis::MatrixEntry>::MPI_DCCols
			(pastis::KmerOps::generate_S(klength, subk_count, alph, parops, tp,
										 local_kmers));
		tp->times["end_main:genS()"] = std::chrono::system_clock::now();
		tu.print_str("Matrix S: ");
		tu.print_str("\nLoad imbalance: " +
					 std::to_string(S->LoadImbalance()) + "\n");
		S->PrintInfo();

		// AxS
		tp->times["start_main:AxS()"] = tp->times["end_main:genS()"];
		proc_log_stream << "INFO: Rank: " << parops->world_proc_rank
						<< " starting AS" << std::endl;
		A = Mult_AnXBn_DoubleBuff<SubKmerIntersectSR_t,
								  pastis::MatrixEntry,
								  PSpMat<pastis::MatrixEntry>::DCCols>(A, *S);
		proc_log_stream << "INFO: Rank: " << parops->world_proc_rank
						<< " done AS" << std::endl;
		tu.print_str("Matrix AS: ");
		tu.print_str("\nLoad imbalance: " +
					 std::to_string(A.LoadImbalance()) + "\n");
		A.PrintInfo();
		tp->times["end_main:AxS()"] = std::chrono::system_clock::now();
		delete S;
	}

	// Output matrix
	tp->times["start_main:(AS)At()"] = std::chrono::system_clock::now();
	proc_log_stream << "INFO: Rank: " << parops->world_proc_rank
					<< " starting AAt" << std::endl;
	PSpMat<pastis::CommonKmers>::MPI_DCCols C =
		Mult_AnXBn_DoubleBuff<KmerIntersectSR_t,pastis::CommonKmers,
							  PSpMat<pastis::CommonKmers>::DCCols>(A, At);
	proc_log_stream << "INFO: Rank: " << parops->world_proc_rank
					<< " done AAt" << std::endl;
	tu.print_str(
		"Matrix AAt: Overlaps after k-mer finding (nnz(C) - diagonal): "
		+ std::to_string(C.getnnz() - seq_count)
		+ "\nLoad imbalance: " + std::to_string(C.LoadImbalance()) + "\n");
	tp->times["end_main:(AS)At()"] = std::chrono::system_clock::now();
	tu.print_str("Matrix B, i.e AAt or ASAt: ");
	C.PrintInfo();

	// the output matrix may be non-symmetric when substitutes are used
	if (add_substitue_kmers)
	{
		auto CT = C;
		CT.Transpose();
		CT.Apply([] (pastis::CommonKmers &arg)
				 {
					 std::swap(arg.first.first, arg.first.second);
					 std::swap(arg.second.first, arg.second.second);
					 return arg;
				 });
		C += CT;

		tu.print_str("Matrix B, i.e AAt or ASAt: (after sym) ");
		C.PrintInfo();
	}
	tu.print_str("Matrix C: ");
	tu.print_str("\nLoad imbalance: " +
				 std::to_string(C.LoadImbalance()) + "\n");

	// -------------------------------------------------------------------------

	// Wait until data distribution is complete
	tp->times["start_main:dfd->wait()"] = std::chrono::system_clock::now();
	if (!dfd->is_ready())
	  	dfd->wait();
	tp->times["end_main:dfd->wait()"] = std::chrono::system_clock::now();

	uint64_t	n_rows			 = dfd->global_count();
	uint64_t	n_cols			 = dfd->global_count();
	int			gr_rows			 = parops->grid->GetGridRows();
  	int			gr_cols			 = parops->grid->GetGridCols();
	int			gr_col_idx		 = parops->grid->GetRankInProcRow();
  	int			gr_row_idx		 = parops->grid->GetRankInProcCol();
	uint64_t	avg_rows_in_grid = n_rows / gr_rows;
	uint64_t	avg_cols_in_grid = n_cols / gr_cols;
	uint64_t	row_offset		 = gr_row_idx * avg_rows_in_grid;
	uint64_t	col_offset		 = gr_col_idx * avg_cols_in_grid;

	// Start pairwise runner
	DistributedPairwiseRunner dpr(dfd, C.seqptr(), &C, afreq,
								  row_offset, col_offset, parops);
	if (!no_align)
	{
		tp->times["start_main:dpr->align()"] = std::chrono::system_clock::now();
		seqan::Blosum62 blosum62(gap_ext, gap_open);
		// align_file += "_Rank_" +
		// 	std::to_string(parops->world_proc_rank) + ".txt";
		// TODO: SeqAn can't work with affine gaps for seed extension
		seqan::Blosum62 blosum62_simple(gap_open, gap_open);
		PairwiseFunction* pf = nullptr;
		uint64_t local_alignments = 0;
		if (xdrop_align)
		{
			pf = new SeedExtendXdrop (blosum62, blosum62_simple,
									  klength, xdrop, seed_count);
			dpr.run_batch(pf, align_file, proc_log_stream, log_freq,
						  ckthr, mosthr * klength, tu);
			local_alignments = static_cast<SeedExtendXdrop*>(pf)->nalignments;
		}
		else if (full_align)
		{
			pf = new FullAligner(blosum62, blosum62_simple);
			dpr.run_batch(pf, align_file, proc_log_stream, log_freq,
						  ckthr, mosthr * klength, tu);
			local_alignments = static_cast<FullAligner*>(pf)->nalignments;
		}
		else if(banded_align)
		{
			pf = new BandedAligner (blosum62, banded_half_width);
			dpr.run_batch(pf, align_file, proc_log_stream, log_freq,
						  ckthr, mosthr * klength, tu);
			local_alignments = static_cast<BandedAligner*>(pf)->nalignments;
		}

		tp->times["end_main:dpr->align()"] = std::chrono::system_clock::now();
		delete pf;
		uint64_t total_alignments = 0;
		MPI_Reduce(&local_alignments, &total_alignments, 1,
				   MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
		tu.print_str("total #alignments: " +
					 std::to_string(total_alignments) + "\n");
	}

	tp->times["start_main:dpr->write_overlaps()"] =
		std::chrono::system_clock::now();
	if (write_overlaps)
	  	dpr.write_overlaps(overlap_file.c_str());
	tp->times["end_main:dpr->write_overlaps()"] =
		std::chrono::system_clock::now();

	tp->times["end_main"] = std::chrono::system_clock::now();

	std::time_t end_prog_time = std::chrono::system_clock::to_time_t(
	  tp->times["end_main"]);
	print_str = "INFO: Program ended on ";
	print_str.append(std::ctime(&end_prog_time));
	tu.print_str(print_str);
	tu.print_str(tp->to_string());

	proc_log_stream.close();

	MPI_Barrier(MPI_COMM_WORLD);
	parops->teardown_parallelism();
	
	return 0;
}


int parse_args(int argc, char **argv) {
  cxxopts::Options options("PASTIS",
                           "A fast, distributed protein aligner");

  options.add_options()
    (CMD_OPTION_INPUT, CMD_OPTION_DESCRIPTION_INPUT,
     cxxopts::value<std::string>())
    (CMD_OPTION_INPUT_SEQ_COUNT, CMD_OPTION_DESCRIPTION_INPUT_SEQ_COUNT,
     cxxopts::value<int>())
    (CMD_OPTION_INPUT_OVERLAP, CMD_OPTION_DESCRIPTION_INPUT_OVERLAP,
     cxxopts::value<uint64_t>())
    (CMD_OPTION_SEED_COUNT, CMD_OPTION_DESCRIPTION_SEED_COUNT,
     cxxopts::value<int>())
    (CMD_OPTION_GAP_OPEN, CMD_OPTION_DESCRIPTION_GAP_OPEN,
     cxxopts::value<int>())
    (CMD_OPTION_GAP_EXT, CMD_OPTION_DESCRIPTION_GAP_EXT,
     cxxopts::value<int>())
    (CMD_OPTION_KMER_LENGTH, CMD_OPTION_DESCRIPTION_KMER_LENGTH,
     cxxopts::value<int>())
    (CMD_OPTION_KMER_STRIDE, CMD_OPTION_DESCRIPTION_KMER_STRID,
     cxxopts::value<int>())
    (CMD_OPTION_OVERLAP_FILE, CMD_OPTION_DESCRIPTION_OVERLAP_FILE,
     cxxopts::value<std::string>())
    (CMD_OPTION_ALIGN_FILE, CMD_OPTION_DESCRIPTION_ALIGN_FILE,
     cxxopts::value<std::string>())
    (CMD_OPTION_NO_ALIGN, CMD_OPTION_DESCRIPTION_NO_ALIGN)
    (CMD_OPTION_FULL_ALIGN, CMD_OPTION_DESCRIPTION_FULL_ALIGN)
    (CMD_OPTION_XDROP_ALIGN, CMD_OPTION_DESCRIPTION_XDROP_ALIGN,
     cxxopts::value<int>())
    (CMD_OPTION_BANDED_ALIGN, CMD_OPTION_DESCRIPTION_BANDED_ALIGN,
     cxxopts::value<int>())
    (CMD_OPTION_IDX_MAP, CMD_OPTION_DESCRIPTION_IDX_MAP,
     cxxopts::value<std::string>())
    (CMD_OPTION_ALPH, CMD_OPTION_DESCRIPTION_ALPH,
     cxxopts::value<std::string>())
    (CMD_OPTION_JOB_NAME_PREFIX, CMD_OPTION_DESCRIPTION_JOB_NAME_PREFIX,
     cxxopts::value<std::string>())
    (CMD_OPTION_SUBS, CMD_OPTION_DESCRIPTION_SUBS,
     cxxopts::value<int>())
    (CMD_OPTION_LOG_FREQ, CMD_OPTION_DESCRIPTION_LOG_FREQ,
     cxxopts::value<int>())
    (CMD_OPTION_AF_FREQ, CMD_OPTION_DESCRIPTION_AF_FREQ,
     cxxopts::value<int>())
	(CMD_OPTION_CKTHR, CMD_OPTION_DESCRIPTION_CKTHR,
     cxxopts::value<int>())
	(CMD_OPTION_MOSTHR, CMD_OPTION_DESCRIPTION_MOSTHR,
     cxxopts::value<float>());

  auto result = options.parse(argc, argv);

  bool is_world_rank0 = parops->world_proc_rank == 0;
  if (result.count(CMD_OPTION_INPUT)) {
    input_file = result[CMD_OPTION_INPUT].as<std::string>();
  } else {
    if (is_world_rank0) {
      std::cout << "ERROR: Input file not specified" << std::endl;
    }
    return -1;
  }

    if (result.count(CMD_OPTION_IDX_MAP)) {
        idx_map_file = result[CMD_OPTION_IDX_MAP].as<std::string>();
    } else {
        if (is_world_rank0) {
            std::cout << "ERROR: Index map file not specified" << std::endl;
        }
        return -1;
    }

  // TODO - fix option parsing
  if (result.count(CMD_OPTION_INPUT_SEQ_COUNT)) {
    seq_count = result[CMD_OPTION_INPUT_SEQ_COUNT].as<int>();
  } else {
    if (is_world_rank0) {
      std::cout << "ERROR: Input sequence count not specified" << std::endl;
    }
    return -1;
  }

  if (result.count(CMD_OPTION_INPUT_OVERLAP)) {
    input_overlap = result[CMD_OPTION_INPUT_OVERLAP].as<uint64_t>();
  } else {
    input_overlap = 10000;
  }

  if (result.count(CMD_OPTION_SEED_COUNT)) {
    seed_count = result[CMD_OPTION_SEED_COUNT].as<int>();
  } else {
    seed_count = 2;
  }

  if (result.count(CMD_OPTION_GAP_OPEN)) {
    gap_open = result[CMD_OPTION_GAP_OPEN].as<int>();
  } else {
    gap_open = -11;
  }

  if (result.count(CMD_OPTION_GAP_EXT)) {
    gap_ext = result[CMD_OPTION_GAP_EXT].as<int>();
  } else {
    gap_ext = -2;
  }

  if (result.count(CMD_OPTION_KMER_LENGTH)) {
    klength = result[CMD_OPTION_KMER_LENGTH].as<int>();
  } else {
    if (is_world_rank0) {
      std::cout << "ERROR: Kmer length not specified" << std::endl;
    }
    return -1;
  }

  if (result.count(CMD_OPTION_KMER_STRIDE)) {
    kstride = result[CMD_OPTION_KMER_STRIDE].as<int>();
  } else {
    if (is_world_rank0) {
      kstride = 1;
    }
  }

  if (result.count(CMD_OPTION_OVERLAP_FILE)) {
    write_overlaps = true;
    overlap_file = result[CMD_OPTION_OVERLAP_FILE].as<std::string>();
  }

  if (result.count(CMD_OPTION_ALIGN_FILE)) {
    align_file = result[CMD_OPTION_ALIGN_FILE].as<std::string>();
  }

  if (result.count(CMD_OPTION_NO_ALIGN)) {
    no_align = true;
  }

  if (result.count(CMD_OPTION_FULL_ALIGN)) {
    full_align = true;
  }

  if (result.count(CMD_OPTION_BANDED_ALIGN)) {
    banded_align = true;
    banded_half_width = result[CMD_OPTION_BANDED_ALIGN].as<int>();
  }

  if (result.count(CMD_OPTION_XDROP_ALIGN)) {
    xdrop_align = true;
    xdrop = result[CMD_OPTION_XDROP_ALIGN].as<int>();
  }

  if (result.count(CMD_OPTION_ALPH)) {
    std::string tmp = result[CMD_OPTION_ALPH].as<std::string>();
    if (tmp == "protein"){
      alph_t = Alphabet::PROTEIN;
    } else {
      if (is_world_rank0) {
        std::cout << "ERROR: Unsupported alphabet type " << tmp << std::endl;
      }
    }
  } else {
    alph_t = Alphabet::PROTEIN;
  }

  if (result.count(CMD_OPTION_SUBS)) {
    add_substitue_kmers = true;
    subk_count = result[CMD_OPTION_SUBS].as<int>();
  }

  // boost::uuids::random_generator gen;
  // boost::uuids::uuid id = gen();
  // if (result.count(CMD_OPTION_JOB_NAME_PREFIX)) {
  //   std::string tmp = result[CMD_OPTION_JOB_NAME_PREFIX].as<std::string>();
  //   job_name = tmp + "_" +  boost::uuids::to_string(id);
  // } else {
  //   job_name = boost::uuids::to_string(id);
  // }

  if (result.count(CMD_OPTION_LOG_FREQ)) {
    log_freq = result[CMD_OPTION_LOG_FREQ].as<int>();
  }

  if (result.count(CMD_OPTION_AF_FREQ)) {
    afreq = result[CMD_OPTION_AF_FREQ].as<int>();
  }

  if (result.count(CMD_OPTION_CKTHR)) {
     ckthr = result[CMD_OPTION_CKTHR].as<int>();
  }

  if (result.count(CMD_OPTION_MOSTHR)) {
     mosthr = result[CMD_OPTION_MOSTHR].as<float>();
  }

  return 0;
}

auto bool_to_str = [](bool b){
  return b ? "True" : "False";
};

void pretty_print_config(std::string &append_to) {
  std::vector<std::string> params = {
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
    "No align (--na)",
    "Full align (--fa)",
    "Xdrop align (--xa)",
    "Banded align (--ba)",
    "Index map (--idxmap)",
    "Alphabet (--alph)",
    "Use substitute kmers (--subs)",
	"Common k-mer threshold (--ckthr)",
	"Max overlap score threshold (--mosthr)"
  };



  std::vector<std::string> vals = {
    input_file,
    std::to_string(seq_count),
    std::to_string(klength),
    std::to_string(kstride),
    std::to_string(input_overlap),
    std::to_string(seed_count),
    std::to_string(gap_open),
    std::to_string(gap_ext),
    !overlap_file.empty() ? overlap_file : "None",
    !align_file.empty() ? align_file : "None",
    !align_file.empty() ? std::to_string(afreq) : "None",
    bool_to_str(no_align),
    bool_to_str(full_align),
    bool_to_str(xdrop_align) + (xdrop_align ? "| xdrop: " + std::to_string(xdrop) : ""),
    bool_to_str(banded_align) + (banded_align ? " | half band: " + std::to_string(banded_half_width) : ""),
    !idx_map_file.empty() ? idx_map_file : "None",
    std::to_string(alph_t),
    bool_to_str(add_substitue_kmers) + (add_substitue_kmers ? " | sub kmers: " + std::to_string(subk_count) : ""),
	std::to_string(ckthr),
	std::to_string(mosthr),
  };

  ushort max_length = 0;
  for (const auto &param : params) {
    if (param.size() > max_length) {
      max_length = static_cast<ushort>(param.size());
    }
  }

  std::string prefix = "  ";
  append_to.append("Parameters...\n");
  for (ushort i = 0; i < params.size(); ++i) {
    std::string param = params[i];
    append_to.append(prefix).append(param).append(": ")
      .append(get_padding(
        static_cast<ushort>(max_length - param.size()), ""))
      .append(vals[i]).append("\n");

//    append_to.append(get_padding(
//        static_cast<ushort>(max_length + 1 - param.row_size()), prefix))
//        .append(param).append(": ")
//        .append(vals[i]).append("\n");
  }
}

std::string get_padding(ushort count, std::string prefix) {
  std::string pad = std::move(prefix);
  for (ushort i = 0; i < count; ++i) {
    pad += " ";
  }
  return pad;
}
