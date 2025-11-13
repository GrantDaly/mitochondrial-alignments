cmake_minimum_required (VERSION 3.28)
project (mitochondrial-alignments CXX)
set(CMAKE_CXX_STANDARD 23)

# make cache variables for install destinations
include(GNUInstallDirs)


#################
# seqan3
################
include(FetchContent)
FetchContent_Declare(seqan3
  GIT_REPOSITORY https://github.com/seqan/seqan3
  GIT_TAG 6b681fb2eae5ab2997d293e99fc6a7f869a20316 # version 3.3.0. Am not trying to upgrade to 3.4.0 if possible, because they deprecated the argument parser and put it as a separate project. I already planned to replace it with getopt. 
  )
FetchContent_MakeAvailable(seqan3)

#################
# hts-wrapper
# Setting up FetchContent_Declare to fetch the dependency. Plan is to override so I can edit both, but it will default to GitHub for a normal installation
################
FetchContent_Declare(HtsWrapper
  GIT_REPOSITORY https://github.com/GrantDaly/hts-wrapper
  )
FetchContent_MakeAvailable(HtsWrapper)

find_package(GSL REQUIRED)

################
# Executables
################

add_executable (inserts-and-cov source/inserts-and-cov.cpp)
target_link_libraries (inserts-and-cov seqan3::seqan3 hts GSL::gsl)

add_executable(identify-numt source/identify-numt.cpp)
target_link_libraries (identify-numt seqan3::seqan3 hts-wrapper hts GSL::gsl)

add_executable(call-variants source/call-variants.cpp source/samplePileup.hpp source/sampleBasePileup.hpp)
target_link_libraries (call-variants seqan3::seqan3 hts GSL::gsl)

add_executable(simulated-reads source/simulated-reads.cpp)
target_link_libraries (simulated-reads seqan3::seqan3 hts-wrapper hts GSL::gsl)

add_executable(mito-align source/mito-align.cpp)
target_link_libraries (mito-align hts-wrapper hts GSL::gsl)
