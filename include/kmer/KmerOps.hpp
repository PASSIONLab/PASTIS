// Created by Saliya Ekanayake on 10/15/19.

#ifndef PASTIS_KMEROPS_HPP
#define PASTIS_KMEROPS_HPP

#include <unordered_set>
#include "MatrixEntry.hpp"
#include "../Types.hpp"
#include "../ParallelOps.hpp"
#include "../Alphabet.hpp"
#include "../ScoreMat.hpp"
#include "../Utils.hpp"
#include "../DistributedFastaData.hpp"
#include "Kmer.hpp"

namespace pastis {
  class KmerOps {
  public:
    static PSpMat<MatrixEntry>::MPI_DCCols generate_A(uint64_t seq_count,
        std::shared_ptr<DistributedFastaData>& dfd, ushort k, ushort s,
        Alphabet &alph, const std::shared_ptr<ParallelOps>& parops,
        const std::shared_ptr<TimePod>& tp, std::unordered_set<Kmer, Kmer>& local_kmers);

    static uint64_t add_kmers(const char *seq, ushort len, uint64_t start_offset,
                              uint64_t end_offset_inclusive, ushort k, ushort s,
                              Alphabet &alp, uvec_64 &lcol_ids, std::vector<MatrixEntry> &lvals,
                              const std::shared_ptr<ParallelOps>& parops,
                              std::unordered_set<Kmer, Kmer>& local_kmers) {

	  pastis::Blosum62 bsm62;

      auto num_kmers = static_cast<ushort>((floor((len - k) * 1.0 / s)) + 1);
      // TODO: Saliya - this can be improved using bit operators
      ushort base = alp.size;
      uint64_t kcode = 0;
      uint64_t count = 0;
      char cap_c;
      std::unordered_set<uint64_t> kmers_in_sequence;
      for (uint64_t i = start_offset;
           i <= ((end_offset_inclusive - k) + 1); i += s) {
        kcode = 0;
        if (count == num_kmers) break;
        std::string kmer_str;
        for (uint64_t j = i; j < i + k; ++j) {
          /*! Efficient than using pow() */
          cap_c = *(seq + j);
          if (cap_c > 96 && cap_c < 123) {
            // small case character, so make it uppercase.
            cap_c = cap_c - 32;
          }
          // kmer_str += cap_c;
		  kmer_str += alp.al_map[cap_c];
          kcode = kcode * base + alp.char_to_code[cap_c];
        }

        ++count;
        lcol_ids.push_back(kcode);
        local_kmers.emplace(kmer_str, kcode, alp, false);

		// cout << kmer_str << " " << kcode << "\n";

        /*! Offset is relative to the sequence start, so unsigned short is
         * good enough. */
		// @OGUZ-EDIT make first entry self score
        lvals.emplace_back(bsm62.self_score(kmer_str),
						   static_cast<ushort &&>(i - start_offset));
      }

      if (count != num_kmers) {
        fprintf(stderr,
                "ERROR: kmerop: add_kmers(): rank: %d, count:%d numk: %d len: %d k: %d s: %d soff: %llu eoffinc: %llu\n",
                parops->world_proc_rank, count, num_kmers, len, k, s,
                start_offset, end_offset_inclusive);
        fflush(stderr);
      }
      return count;
    }

    static PSpMat<MatrixEntry>::MPI_DCCols generate_S(
        ushort k, uint64_t subk_count,
        Alphabet &alph, std::shared_ptr<ParallelOps> &parops,
        std::shared_ptr<TimePod> &tp,
        std::unordered_set<Kmer, Kmer>& local_kmers);

    static void generate_S_ineff(uvec_64& lcol_ids, ushort k, Alphabet& alph){


    }


    /*! This is buggy */
//    static void generate_S(
//        ushort k, Alphabet& alph, const std::shared_ptr<ParallelOps>& parops)
//    {
//      auto alph_size = alph.size;
//      auto world_rank = parops->world_proc_rank;
//      auto world_size = parops->world_procs_count;
//
//      uint64_t uni_size = std::pow((uint64_t)alph_size, k);
//      uint64_t q = uni_size / world_size;
//      uint64_t r = uni_size - (q * world_size);
//
//      uint64_t kmer_start_offset = q * world_rank + (world_rank < r ? world_rank : r);
//      uint64_t kmer_count = world_rank < r ? q+1 : q;
//
//      std::stringstream ss;
//      ss << "Rank: " << world_rank << " kso:" << kmer_start_offset << " kcount:" << kmer_count << std::endl;
//      std::cout << ss.str();
//
//      std::vector<ushort> sids(k);
//      uint64_t tmp = kmer_start_offset;
//      std::string kmer_str;
//      for (size_t i = 0; i < k; ++i){
//        q = tmp / alph_size;
//        r = tmp - (q*alph_size);
//        kmer_str = alph[r] + kmer_str;
//        sids[i] = r;
//        tmp = q;
//      }
//
//      Kmer root(kmer_str, alph);
//      ss.str(std::string());
//      ss << world_rank << ":" << root;
//      std::cout << ss.str();
//
//      std::queue<Kmer> qew;
//      qew.push(root);
//      std::vector<Kmer> nbrs;
//      while (!qew.empty() && nbrs.size() < kmer_count) {
//        Kmer kmer = qew.front();
//        qew.pop();
//        nbrs.push_back(kmer);
//        for (auto &fid : kmer.get_free_idxs()) {
//          int count = 0;
//          while ((nbrs.size() + qew.size()) < kmer_count &&
//                 (++count + sids[fid]) < alph_size) {
//            Kmer next = kmer.substitute(fid, alph[sids[fid] + count], 0, alph);
//            qew.push(next);
//          }
//        }
//      }
//
//      ss.str(std::string());
//      ss << "Rank: " << world_rank << " " << nbrs.size() << std::endl;
//      for (auto& kmer : nbrs){
//        ss << "  " << kmer;
//      }
//
//      std::cout << ss.str();
//    }

  };
}


#endif //PASTIS_KMEROPS_HPP
