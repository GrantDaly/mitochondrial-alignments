cmake_minimum_required (VERSION 3.4)
project (mitochondrial_processing CXX)
set(CMAKE_CXX_STANDARD 23)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# make cache variables for install destinations
include(GNUInstallDirs)



# add seqan3 to search path
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../seqan3/build_system")
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../hts-wrapper/build/cmake")
 
# require seqan3 with a version between >=3.0.0 and <4.0.0
find_package (seqan3 3.0 REQUIRED)
# supress warning about hardware destructive interference size. Taken from seqan Github. 
#As far as I can tell this could only impact ABI compatibility and this software is not inteneded to be portable in the compiled state.
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        add_compile_options ( -Wall -Wextra -Wpedantic)
        # GCC12 and above: Disable warning about std::hardware_destructive_interference_size not being ABI-stable.
        if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 12)
            add_compile_options ("-Wno-interference-size")
        endif ()
    endif ()


find_package(hts-wrapper REQUIRED)
find_package(GSL REQUIRED)
#find_package(Boost)

# build app with seqan3


#add_executable (inserts-and-cov inserts-and-cov.cpp)
#target_link_libraries (inserts-and-cov seqan3::seqan3 hts gsl)

add_executable (inserts-and-cov inserts-and-cov.cpp)
target_link_libraries (inserts-and-cov seqan3::seqan3 hts GSL::gsl)

add_executable(identify-numt identify-numt.cpp)
target_link_libraries (identify-numt seqan3::seqan3 hts-wrapper hts GSL::gsl)

add_executable(call-variants call-variants.cpp samplePileup.hpp sampleBasePileup.hpp)
target_link_libraries (call-variants seqan3::seqan3 hts GSL::gsl)

add_executable(simulated-reads simulated-reads.cpp)
target_link_libraries (simulated-reads seqan3::seqan3 hts-wrapper hts GSL::gsl)

add_executable(mito-align mito-align.cpp)
target_link_libraries (mito-align hts-wrapper hts GSL::gsl)
