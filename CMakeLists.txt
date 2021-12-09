cmake_minimum_required(VERSION 3.11)
set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")
project(pastis LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 14)

include_directories(
	${PROJECT_SOURCE_DIR}/seqan/include 
	${PROJECT_SOURCE_DIR}/CombBLAS/install/include
	${PROJECT_SOURCE_DIR}/CombBLAS/Applications
	${PROJECT_SOURCE_DIR}/ADEPT
	$ENV{ROCM_PATH}/include
	)

link_directories(
    ${PROJECT_SOURCE_DIR}/CombBLAS/install/lib
    ${PROJECT_SOURCE_DIR}/ADEPT/build
    $ENV{ROCM_PATH}/lib
    )

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DSEQAN_ENABLE_DEBUG=0")
add_compile_definitions(NDEBUG)

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

FILE(GLOB PASTIS_CXX_SRCS src/*.cpp src/kmer/*.cpp src/align/*.cpp)
FILE(GLOB PASTIS_CXX_INCS inc/*.hpp inc/kmer/*.hpp inc/align/*.hpp)
add_executable(pastis
	${PASTIS_CXX_SRCS}
	${PASTIS_CXX_INCS}
	)

# target_link_libraries(pastis PRIVATE -fsanitize=address)
target_link_libraries(pastis PUBLIC CombBLAS GraphGenlib Usortlib)
if (MPI_FOUND)
	message(STATUS "MPI found.")
    include_directories(${MPI_INCLUDE_PATH})
    target_link_libraries(pastis PUBLIC ${MPI_LIBRARIES})
endif()

# ADEPT
target_link_libraries(pastis PUBLIC adept)

if(OPENMP_FOUND OR OpenMP_CXX_FOUND)
    # Set THREADED if OpenMP is found
    target_compile_definitions(pastis PUBLIC THREADED)
    target_compile_options(pastis PUBLIC "${OpenMP_CXX_FLAGS}")
    target_link_libraries(pastis PUBLIC "${OpenMP_CXX_FLAGS}")
endif()

target_link_libraries(pastis PUBLIC amdhip64 hsa-runtime64)

