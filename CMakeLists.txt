cmake_minimum_required(VERSION 3.10)
project(ConsMLP_lw VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories(${PROJECT_SOURCE_DIR}/include)

if(MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /DNDEBUG /GL /fp:fast")
    set(CMAKE_CXX_FLAGS_DEBUG   "/Zi /Od /W3")
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
else()
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native -flto -funroll-loops -ffast-math")
    set(CMAKE_CXX_FLAGS_DEBUG   "-g -O0 -Wall -Wextra")
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
    message(STATUS "Build type not specified: Using Release as default")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

find_package(Threads REQUIRED)

set(LIB_SOURCES
    src/datastructures/Hypergraph.cpp
    src/datastructures/HypergraphHierarchy.cpp
    src/partitioning/Partition.cpp
    src/partitioning/PartitionConstraints.cpp
    src/partitioning/PartitionMetrics.cpp
    src/partitioning/Partitioner.cpp
    src/partitioning/MultilevelPartitioner.cpp
    src/coarsening/Coarsener.cpp
    src/coarsening/MultilevelCoarsener.cpp
    src/refinement/Refiner.cpp
    src/utils/HgrParser.cpp
)

add_executable(ConsMLP_lw src/ConsMLP_lw.cpp ${LIB_SOURCES})
target_link_libraries(ConsMLP_lw PRIVATE Threads::Threads)
