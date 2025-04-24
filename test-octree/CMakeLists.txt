cmake_minimum_required(VERSION 3.20)
project(test-octree VERSION 1.0.0 LANGUAGES C)

set(AMRACUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../build/install/)

# Find MPI package
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# compiler flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -march=native -fopenmp -pedantic -Wall -Wextra -Wshadow -fdiagnostics-color=always")




set(TARGET_NAME "test-octree")


add_executable(${TARGET_NAME} main.c)

target_link_libraries(${TARGET_NAME} PRIVATE MPI::MPI_C)


# Add the include directory to the search path
target_include_directories(${TARGET_NAME} PRIVATE ${AMRACUT_PATH}/include)
target_link_libraries(${TARGET_NAME} PRIVATE ${AMRACUT_PATH}/lib/libamracut.so)


