cmake_minimum_required(VERSION 3.20)
project(amracut VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)


if(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE Release)

endif()



# compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -fopenmp -pedantic -Wall -Wextra -Wshadow -fdiagnostics-color=always")

set(CMAKE_CXX_FLAGS_DEBUG "-g -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(TARGET_NAME "amracut")
add_library(${TARGET_NAME} SHARED)

if(NOT DEFINED AMRACUT_INTEGER_WIDTH)
    set(AMRACUT_INTEGER_WIDTH 32)
endif()

if(AMRACUT_INTEGER_WIDTH EQUAL 64)
    set(AMRACUT_UINT_T "uint64_t")
    set(AMRACUT_SINT_T "int64_t")
    set(AMRACUT_OCTREE_NO_ELEM "UINT64_MAX")
elseif(AMRACUT_INTEGER_WIDTH EQUAL 32)
    set(AMRACUT_UINT_T "uint32_t")
    set(AMRACUT_SINT_T "int32_t")
    set(AMRACUT_OCTREE_NO_ELEM "UINT32_MAX")
elseif(AMRACUT_INTEGER_WIDTH EQUAL 16)
    set(AMRACUT_UINT_T "uint16_t")
    set(AMRACUT_SINT_T "int16_t")
    set(AMRACUT_OCTREE_NO_ELEM "UINT16_MAX")
else()
    message(FATAL_ERROR "Unsupported AMRACUT_INTEGER_WIDTH: ${AMRACUT_INTEGER_WIDTH}. Allowed values: 64, 32, 16")
endif()


message("Using AMRACUT_INTEGER_WIDTH = ${AMRACUT_INTEGER_WIDTH}")

target_compile_definitions(${TARGET_NAME} PUBLIC AMRACUT_INTEGER_WIDTH=${AMRACUT_INTEGER_WIDTH})



# Find MPI package
# include(FindMPI)
# target_include_directories(${TARGET_NAME} PRIVATE ${MPI_INCLUDE_PATH})

find_package(MPI REQUIRED)
target_link_libraries(${TARGET_NAME} PUBLIC MPI::MPI_CXX)



add_subdirectory(include)

add_subdirectory(src)

# Add the include directory to the search path
target_include_directories(${TARGET_NAME} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/include)
target_include_directories(${TARGET_NAME} PRIVATE src)


# Set default install directory to build_dir/install if not provided
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "Installation Directory" FORCE)
endif()


# Install the shared library
install(TARGETS ${TARGET_NAME}
    LIBRARY DESTINATION lib
)




# Installing public headers

# Use file(GLOB) to get all .h files in the 'data' directory
file(GLOB gen_header_files "${CMAKE_BINARY_DIR}/include/*.h")

# Install the files that match the pattern
install(FILES ${gen_header_files} DESTINATION include)


