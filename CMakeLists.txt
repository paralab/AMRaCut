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

target_compile_definitions(${TARGET_NAME} PUBLIC AMRACUT_INTEGER_WIDTH=${AMRACUT_INTEGER_WIDTH})


# add_compile_definitions(AMRACUT_INTEGER_WIDTH=${AMRACUT_INTEGER_WIDTH})

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

# target_include_directories(${TARGET_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/src/dist-graph)


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

# set(DRIVER_TARGET_NAME "test")


# add_executable(${DRIVER_TARGET_NAME} test/main.c)

# # target_link_libraries(${DRIVER_TARGET_NAME} PRIVATE MPI::MPI_C)


# # Add the include directory to the search path
# target_include_directories(${DRIVER_TARGET_NAME} PRIVATE ${CMAKE_BINARY_DIR}/include)
# target_link_libraries(${DRIVER_TARGET_NAME} PRIVATE ${TARGET_NAME})
# target_link_libraries(${DRIVER_TARGET_NAME} PUBLIC MPI::MPI_C)
