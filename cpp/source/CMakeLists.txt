cmake_minimum_required (VERSION 3.4)
project (mitochondrial_processing CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# make cache variables for install destinations
include(GNUInstallDirs)

# add seqan3 to search path
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../seqan3/build_system")
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../hts-wrapper/build/cmake")
 
# require seqan3 with a version between >=3.0.0 and <4.0.0
find_package (seqan3 3.0 REQUIRED)

find_package(hts-wrapper REQUIRED)
find_package(Boost)

# build app with seqan3


#add_executable (inserts-and-cov inserts-and-cov.cpp)
#target_link_libraries (inserts-and-cov seqan3::seqan3 hts gsl)

add_executable (inserts-and-cov inserts-and-cov.cpp)
target_link_libraries (inserts-and-cov seqan3::seqan3 hts gsl)

add_executable(identify-numt identify-numt.cpp)
target_link_libraries (identify-numt seqan3::seqan3 hts gsl)

add_executable(call-variants call-variants.cpp samplePileup.hpp sampleBasePileup.hpp)
target_link_libraries (call-variants seqan3::seqan3 hts gsl)
