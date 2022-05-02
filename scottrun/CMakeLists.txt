cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(msuboltz_project)
set(GITHOME_MSU $ENV{GITHOME_MSU})

add_custom_target(extern_msu_commonutils COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/commonutils/software)
add_custom_target(extern_msu_sampler COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/msu_sampler/software)
add_custom_target(extern_msu_boltzmann COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/msu_boltzmann/software)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -std=c++11 -mfpmath=sse")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL 1.15 REQUIRED)
find_package(Eigen3 REQUIRED)

include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${GITHOME_MSU}/commonutils/software/include  ${GITHOME_MSU}/msu_sampler/software/include  ${GITHOME_MSU}/msu_boltzmann/software/include)

link_directories(${GITHOME_MSU}/msu_sampler/software/lib ${GITHOME_MSU}/msu_boltzmann/software/lib ${GITHOME_MSU}/commonutils/software/lib)

find_library(${GITHOME_MSU}/commonutils/software/lib msu_commonutils)
find_library(${GITHOME_MSU}/msu_sampler/software/lib msu_sampler)
find_library(${GITHOME_MSU}/msu_boltzmann/software/lib msu_boltzmann)

find_library(${GITHOME_MSU}/msu_sampler/software/lib msu_sampler)
find_library(${GITHOME_MSU}/commonutils/software/lib msu_commonutils)
find_library(${GITHOME_MSU}/commonutils/software/lib msu_boltzmann)

set(SOURCE_FILES ./msuboltz.cc)

add_executable(msuboltz ${SOURCE_FILES})

add_dependencies(msuboltz extern_msu_sampler extern_msu_commonutils extern_msu_boltzmann)

target_link_libraries(msuboltz msu_boltzmann msu_commonutils msu_sampler ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY} ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a ${GITHOME_MSU}/msu_boltzmann/software/lib/libmsu_boltzmann.a ${GITHOME_MSU}/msu_sampler/software/lib/libmsu_sampler.a)
