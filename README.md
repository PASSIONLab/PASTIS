
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


2. Build CombBLAS; the following commands can be used to build and install CombBLAS:
```
cd CombBLAS
mkdir build
mkdir install
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install ../
make -j4
make install         
```
## GPU support
PASTIS has GPU support for NVIDIA and AMD GPUs. Switch to CUDA or HIP branch to run PASTIS on systems with accelerators. The alignment option should be provided as ```--absw``` when running on GPUs as it is the alignment option supported for GPUs. To compile the ADEPT library, use:
  * cd ADEPT
  * mkdir build 
  * cd build
  * cmake ..
  * make 

In building the HIP branch, to build ADEPT you need to use the custom cmake file CMakeLists.txt-ADEPT provided in the repo. Copy this file into ADEPT directory by ```cp CMakeLists.txt-ADEPT ADEPT/CMakeLists.txt``` and follow the instructions provided above. In order to use HIP compiler, use ```cmake .. -DCMAKE_CXX_COMPILER=hipcc``` to configure make.
  
## Build and Test PASTIS
To build PASTIS after cloning or downloading the source, use:
  * mkdir build install
  * cd build
  * cmake .. -DCMAKE_INSTALL_PREFIX=../install
  * make install
To run tests, use:
  * ctest -V
This will run PASTIS with the fasta file provided under tests directory, and
produce output alignment information in the same directory.

  
## Run PASTIS

You can run PASTIS in parallel by specifying the number of processes to the mpirun or mpiexec command. The number of processes must be perfect square value.

The main parameters of PASTIS are as follows:
- ```-i <string>```: Input FASTA file.
- ```-c <integer>```: Number of sequences in the FASTA file.
- ```-k <integer>```: K-mer length.
- ```--na```: Do not perform alignment.
- ```--sfa```: Seqan Smith-Waterman alignment.
- ```--sxa <integer>```: Seqan X-drop alignment with the indicated drop value.
- ```--absw```: ADEPT Striped Smith-Waterman (only available for GPUs).
- ```--af <string>```: Final similarity matrix in matrix market file format.

The additional parameters of PASTIS are as follows:
- ```-O <integer>```: Fasta file parallel read overlap in bytes. ```[default: 10000]```
- ```--sc <integer>```: Maximum seed count. ```[default: 2]```
- ```-s <integer>```: K-mer stride. ```[default: 1]```
- ```-g <integer>```: Gap open penalty (negative). ```[default: -11]```
- ```-e <integer>```: Gap extension penalty (negative). ```[default: -2]```
- ```--ckthr <integer>```: Common k-mer threshold. The sequence pairs that have
  less or equal to this number of common k-mers are not aligned. ```[default: 0]```
- ```--mosthr <float>```: Maximum overlap score threshold. The maximum overlap
  score is given by the highest overlap region of k-mers divided by the number
  k-mer length. The pairs whose such score is less than or equal to ```mosthr```
  are not aligned. ```[default: 0]```
- ```--alph <string>```: Reduced alphabet to use. Default is not to use any kind of reduced alphabet. Choose from ```{pdefault, murphy10, dssp10, gbmr10, td10, diamond}```.
- ```--bsz <integer>```: Batch alignment size. ```[default: 1e7]```
- ```--br <integer>```: Block multiplication row dimension. ```[default: 1]```
- ```--bc <integer>```: Block multiplication column dimension. ```[default: 1]```

Example run:
```mpirun -N 4 ./build/pastis -i ./tests/scope_77k.fa -c 77040  --af sim_mat.mtx --ckthr 1 --sxa 49```


## Notes for running PASTIS on large datasets
PASTIS supports a blocking mode for large datasets to perform similarity search in a blocked manner and save memory. In default, this mode is not enabled. If you run into out of memory errors, provide the number of blocks you want PASTIS to use with the options  ```--br``` and ```--bc```, where the target number of blocks is equal to the multiplication of these two parameters.


## Notes for running PASTIS on NERSC Cori
The necessary modules for running PASTIS on Cori are as follows:
* module swap PrgEnv-intel PrgEnv-gnu
* module load cmake

Build PASTIS with the wrappers:
* cmake .. -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC
* make


## Citation

If you use PASTIS in your research, please cite our paper describing the algorithm and the implementation:

  * Oguz Selvitopi*, Saliya Ekanayake*, Giulia Guidi, Georgios Pavlopoulos, Ariful Azad, and Aydın Buluç. Distributed Many-to-Many Protein Sequence Alignment Using Sparse Matrices. Proceedings of the International Conference for High Performance Computing, Networking, Storage, and Analysis (SC’20), 2020. (*:joint first authors) [pdf](https://people.eecs.berkeley.edu/~aydin/PASTIS-SC20.pdf)
