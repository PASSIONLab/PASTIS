
# PASTIS: Distributed Many-to-Many Protein Sequence Alignment Using Sparse Matrices

PASTIS Copyright (c) 2020, The Regents of the University of California,
through Lawrence Berkeley National Laboratory (subject to receipt of
any required approvals from the U.S. Dept. of Energy). All rights reserved.

If you have questions about your rights to use or distribute this software,
please contact Berkeley Lab's Intellectual Property Office at
IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department
of Energy and the U.S. Government consequently retains certain rights.  As
such, the U.S. Government has been granted for itself and others acting on
its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
Software to reproduce, distribute copies to the public, prepare derivative 
works, and perform publicly and display publicly, and to permit others to do so.

## Prerequisites

1. Operating System.
  * PASTIS is tested and known to work on the following operating systems.
    *  SUSE Linux Enterprise Server 15.
    *  Ubuntu 14.10.
    *  MacOS.
    
2. GCC/G++ version 8.2.0 or above.

3. CMake 3.11 or above.

## Dependencies

1. PASTIS requires CombBLAS and SeqAn libraries. These are included within PASTIS as
submodules. To get them all, use:
```
git clone --recursive https://github.com/PASSIONLab/PASTIS
```

If you already happen to have those libraries, they are needed within the PASTIS
directory. So, link the installation directories under names "CombBLAS" and
"seqan". Note that PASTIS uses master branch of CombBLAS and develop branch of
SeqAn.


2. Build CombBLAS:
  * The following commands can be used to build and install CombBLAS:
  ```
    cd CombBLAS
    mkdir build
    mkdir install
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install ../
    make -j4
    make install         
  ```
3. SeqAn is a header only library, so there's no need to build it.
  
## Build and Test PASTIS

To build PASTIS after cloning or downloading the source, use:
  * mkdir build_release
  * cd build_release
  * cmake ..
  * make

To run tests, use:
  * ctest -V

This will run PASTIS with the fasta file provided under tests directory, and
produce output alignment information in the same directory.

  
## Run PASTIS

You can run PASTIS in parallel by specifying the number of processes to the mpirun or mpiexec command. The number of processes must be perfect square value.

The parameters and options of PASTIS are as follows:
- ```-i <string>```: Input FASTA file.
- ```-c <integer>```: Number of sequences in the FASTA file.
- ```--sc <integer>```: Seed count. ```[default: 2]```
- ```-k <integer>```: K-mer length.
- ```-s <integer>```: K-mers stride ```[default: 1]```
- ```--subs <integer>```: Number of substitute K-mers. 
- ```-g <integer>```: Gap open penalty (negative). ```[default: -11]```
- ```-e <integer>```: Gap extension penalty (negative). ```[default: -2]```
- ```-O <integer>```: Number of bytes to overlap when reading the input file in parallel. ```[default: 10000]```
- ```--na```: Do not perform alignment.
- ```--fa```: Smith-Waterman alignment.
- ```--xa <integer>```: X-drop alignment with the indicated drop value.
- ```--ba <integer>```: Banded alignment with the indicated band size.
- ```--af <string>```: Output file to write alignment information. 
- ```--idxmap <string>```: Output file for input sequences to ids used in PASTIS.
- ```--ckthr <integer>```: Common k-mer threshold. The sequence pairs that have
  less or equal to this number of common k-mers are not aligned. ```[default: 0]```
- ```--mosthr <float>```: Maximum overlap score threshold. The maximum overlap
  score is given by the highest overlap region of k-mers divided by the number
  k-mer length. The pairs whose such score is less than or equal to ```mosthr```
  are not aligned. ```[default: 0]```


## Notes for running PASTIS on NERSC Cori

The necessary modules for running PASTIS on Cori are as follows:
* module swap PrgEnv-intel PrgEnv-gnu
* module load cmake

Make sure to pass the correct MPI wrappers to successfully run the tests:
* cmake -DMPIEXEC_EXECUTABLE=/usr/bin/srun ..
* make
* ctest -V


## Citation

If you use PASTIS in your research, please cite our paper describing the algorithm and the implementation:

  * Oguz Selvitopi*, Saliya Ekanayake*, Giulia Guidi, Georgios Pavlopoulos, Ariful Azad, and Aydın Buluç. Distributed Many-to-Many Protein Sequence Alignment Using Sparse Matrices. Proceedings of the International Conference for High Performance Computing, Networking, Storage, and Analysis (SC’20), 2020. (*:joint first authors) [pdf](https://people.eecs.berkeley.edu/~aydin/PASTIS-SC20.pdf)
