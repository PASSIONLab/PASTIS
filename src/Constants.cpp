// Created by Saliya Ekanayake on 12/17/18.


#include "../inc/Constants.hpp"



namespace
pastis
{

const char *CMD_OPTION_INPUT = "i";
const char *CMD_OPTION_DESCRIPTION_INPUT = "The input FASTA file.";

const char *CMD_OPTION_INPUT_SEQ_COUNT = "c";
const char *CMD_OPTION_DESCRIPTION_INPUT_SEQ_COUNT = "The sequence count in the input sequence file.";

const char *CMD_OPTION_INPUT_OVERLAP = "O";
const char *CMD_OPTION_DESCRIPTION_INPUT_OVERLAP = "Number of bytes to overlap when reading the input file in parallel.";

const char *CMD_OPTION_SEED_COUNT = "sc";
const char *CMD_OPTION_DESCRIPTION_SEED_COUNT = "Seed count.";

const char *CMD_OPTION_GAP_OPEN = "g";
const char *CMD_OPTION_DESCRIPTION_GAP_OPEN = "Gap open penalty (negative value).";

const char *CMD_OPTION_GAP_EXT = "e";
const char *CMD_OPTION_DESCRIPTION_GAP_EXT = "Gap extension penalty (negative value).";

const char *CMD_OPTION_KMER_LENGTH = "k";
const char *CMD_OPTION_DESCRIPTION_KMER_LENGTH = "Kmer length.";

const char *CMD_OPTION_KMER_STRIDE = "s";
const char *CMD_OPTION_DESCRIPTION_KMER_STRID = "Kmer stride.";

const char *CMD_OPTION_OVERLAP_FILE = "of";
const char *CMD_OPTION_DESCRIPTION_OVERLAP_FILE = "Overlap file.";

const char *CMD_OPTION_ALIGN_FILE = "af";
const char *CMD_OPTION_DESCRIPTION_ALIGN_FILE = "Alignment file.";

const char *CMD_OPTION_NO_ALIGN = "na";
const char *CMD_OPTION_DESCRIPTION_NO_ALIGN = "Flag to indicate not to performa alignments.";

const char *CMD_OPTION_FULL_ALIGN = "sfa";
const char *CMD_OPTION_DESCRIPTION_FULL_ALIGN = "Flag to indicate full alignment";

const char *CMD_OPTION_XDROP_ALIGN = "sxa";
const char *CMD_OPTION_DESCRIPTION_XDROP_ALIGN = "Flag to indicate seed-and-extend using xdrop alignment";

const char *CMD_OPTION_BANDED_ALIGN = "sba";
const char *CMD_OPTION_DESCRIPTION_BANDED_ALIGN = "Flag to indicate banded alignment";

const char *CMD_OPTION_GPUBSW_ALIGN = "absw";
const char *CMD_OPTION_DESCRIPTION_GPUBSW_ALIGN = "Flag to indicate GPU-based Banded Smith-Waterman";

const char *CMD_OPTION_IDX_MAP = "idxmap";
const char *CMD_OPTION_DESCRIPTION_IDX_MAP = "The file path to write the global sequence indices to original global sequence indices mapping. ";

const char *CMD_OPTION_ALPH = "alph";
const char *CMD_OPTION_DESCRIPTION_ALPH = "The alphabet to use. Valid choices are [pdefault|murphy10|dssp10|gbmr10|td10|diamond].";

const char *CMD_OPTION_SUBS = "subs";
const char *CMD_OPTION_DESCRIPTION_SUBS = "Number of substitute kmers.";

const char *CMD_OPTION_JOB_NAME_PREFIX = "jp";
const char *CMD_OPTION_DESCRIPTION_JOB_NAME_PREFIX = "Job name prefix.";

const char *CMD_OPTION_LOG_FREQ = "lf";
const char *CMD_OPTION_DESCRIPTION_LOG_FREQ = "Log frequency.";

const char *CMD_OPTION_AF_FREQ = "afreq";
const char *CMD_OPTION_DESCRIPTION_AF_FREQ = "Alignment write frequency in number of lines.";

const char *CMD_OPTION_CKTHR = "ckthr";
const char *CMD_OPTION_DESCRIPTION_CKTHR = "Perform alignments that have more than this number of common k-mers.";

const char *CMD_OPTION_MOSTHR = "mosthr";
const char *CMD_OPTION_DESCRIPTION_MOSTHR = "Perform alignments that have more than this overlap score per base.";

const char *CMD_OPTION_ALN_BATCH_SZ = "bsz";
const char *CMD_OPTION_DESCRIPTION_ALN_BATCH_SZ = "Number of pairwise alignments to perform in a batch.";

const char *CMD_OPTION_SIM_BR = "br";
const char *CMD_OPTION_DESCRIPTION_SIM_BR = "Block row dimension.";

const char *CMD_OPTION_SIM_BC = "bc";
const char *CMD_OPTION_DESCRIPTION_SIM_BC = "Block col dimension.";

const char *CMD_OPTION_LB = "lb";
const char *CMD_OPTION_DESCRIPTION_LB = "Load balancing (default: idx) [trg | idx].";

const char *CMD_OPTION_STATS = "stats";
const char *CMD_OPTION_DESCRIPTION_STATS = "Print runtime and various statistics.";

}

