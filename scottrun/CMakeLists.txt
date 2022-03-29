cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
set(GITHOME_MSU $ENV{GITHOME_MSU})

add_custom_target(extern_commonutils COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/commonutils/software)
add_custom_target(extern_msu_sampler COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/msu_sampler/software)
#add_custom_target(extern_msu_boltzmann COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/msu_boltzmann/software)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Wold-style-cast -Werror=switch -std=c++11 -mfpmath=sse")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL 1.15 REQUIRED)
find_package(Eigen3 REQUIRED)

include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${GITHOME_MSU}/commonutils/software/include  ${GITHOME_MSU}/msu_sampler/software/include  ${GITHOME_MSU}/msu_boltzmann/software/include)

link_directories(${GITHOME_MSU}/commonutils/software/lib ${GITHOME_MSU}/msu_sampler/software/lib ${GITHOME_MSU}/msu_boltzmann/software/lib )
message("Libraries: b3d ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY}" )

find_library(${GITHOME_MSU}/commonutils/software/lib commonutils)
find_library(${GITHOME_MSU}/msu_sampler/software/lib msu_samper)
find_library(${GITHOME_MSU}/msu_boltzmann/software/lib msu_boltzmann)

set(SOURCE_FILES ./msuboltz.cc)

add_executable(msuboltz ${SOURCE_FILES})
add_dependencies(msuboltz extern_commonutils extern_msu_sampler)
# extern_msu_boltzmann)

target_link_libraries(msuboltz ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY} ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a ${GITHOME_MSU}/msu_sampler/software/lib/libmsu_sampler.a)
# ${GITHOME_MSU}/software/lib/libmsu_boltzmann.a)
