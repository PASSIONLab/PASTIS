// Created by Saliya Ekanayake on 10/15/19.

#ifndef PASTIS_KMERINTERSECTSR_HPP
#define PASTIS_KMERINTERSECTSR_HPP

#include "../ParallelOps.hpp"
namespace pastis {
  template<typename IN, typename OUT>
  struct KmerIntersect {
    static OUT id() {
      OUT a;
      return a;
    }

    static bool returnedSAID() { return false; }

    static OUT add(const OUT &arg1, const OUT &arg2) {
      OUT res(arg1.count + arg2.count);
      // TODO: perhaps improve this late with something that'll check how far
      // apart are the kmers.
      res.first.first = arg1.first.first;
      res.first.second = arg1.first.second;
      res.second.first = arg2.first.first;
      res.second.second = arg2.first.second;
	  res.score = std::max(arg1.score, arg2.score);
      return res;
    }

//    /* This doesn't make sense */
//    static IN add(const OUT &arg1, const IN &arg2) {
//      IN res(arg2.cost, arg2.offset);
//      return res;
//    }

    static OUT multiply(const IN &arg1, const IN &arg2) {
      OUT a;
      a.first.first = arg1.offset;
      a.first.second = arg2.offset;
	  a.score = std::max(arg1.cost, arg2.cost);
      return a;
    }

    static void axpy(IN a, const IN &x, OUT &y) {
      y = add(y, multiply(a, x));
    }

    static MPI_Op mpi_op() {
      static MPI_Op mpiop;
      static bool exists = false;
      if (exists)
        return mpiop;
      else {
        MPI_Op_create(MPI_func, true, &mpiop);
        exists = true;
        return mpiop;
      }
    }

    static void
    MPI_func(void *invec, void *inoutvec, int *len, MPI_Datatype *datatype) {
      for (int i = 0; i < *len; ++i) {
        *((OUT) inoutvec + i) = add(*((OUT) invec + i), *((OUT) inoutvec + 1));
      }

    }
  };
}
#endif //PASTIS_KMERINTERSECTSR_HPP
