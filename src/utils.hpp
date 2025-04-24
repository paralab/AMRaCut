/**
 * @file    utils.hpp
 * @author  
 * @brief   Utility functions
 * @version 0.1
 * @date 2024-10-25
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef AMRACUT_UTILS_H__
#define AMRACUT_UTILS_H__

#include "amracut.h"
#include <sstream>
#include <iostream>
#include <vector>

// TODO: allow individual integer width customization for distance and label


using graph_indexing_t = amracut_uint_t;
using bfs_distance_t   = amracut_uint_t;
using bfs_label_t      = amracut_uint_t;

#define AMRACUT_BFS_INFINITY std::numeric_limits<amracut_uint_t>::max()
#define AMRACUT_BFS_NO_LABEL std::numeric_limits<amracut_uint_t>::max()



using amracut_real_t = float;



/**
 * @brief << operator overloading for printing a std::pair type
 * 
 * @tparam T first type
 * @tparam U second type
 * @param os output stream
 * @param p pair
 * @return std::ostream& 
 */
template<typename T, typename U>
std::ostream& operator<<(std::ostream& os, const std::pair<T, U>& p) {
    os << "(" << p.first << ", " << p.second << ")";
    return os;
}

namespace amracut
{
  template <typename T>
  class Mpi_datatype;

  /**
   * @brief Template Classes to convert C++ types to MPI data types
   * 
   */
  #define HS_MPIDATATYPE(CTYPE, MPITYPE)	\
    template <>  \
      class Mpi_datatype<CTYPE> \
      {  \
          public: \
        static MPI_Datatype value() {\
                            return MPITYPE;\
                          } \
      };

      HS_MPIDATATYPE(short,          MPI_SHORT)
      HS_MPIDATATYPE(int,            MPI_INT)
      HS_MPIDATATYPE(long,           MPI_LONG)
      HS_MPIDATATYPE(unsigned short, MPI_UNSIGNED_SHORT)
      HS_MPIDATATYPE(unsigned int,   MPI_UNSIGNED)
      HS_MPIDATATYPE(unsigned long,  MPI_UNSIGNED_LONG)
      HS_MPIDATATYPE(float,          MPI_FLOAT)
      HS_MPIDATATYPE(double,         MPI_DOUBLE)
      HS_MPIDATATYPE(long double,    MPI_LONG_DOUBLE)
      HS_MPIDATATYPE(long long,    MPI_LONG_LONG_INT)
      HS_MPIDATATYPE(char,    MPI_CHAR)
      HS_MPIDATATYPE(unsigned char,    MPI_UNSIGNED_CHAR)


  #undef HS_MPIDATATYPE

  /**
   * @brief Template class to convert a std::pair<T1,T2> to MPI data type.
   * 
   * @tparam T1 : Type 1
   * @tparam T2 : Type 2
   */
  template <typename T1, typename T2>
  class Mpi_pairtype
  {
  public:
    static MPI_Datatype value()
    {
      static bool first = true;
      static MPI_Datatype datatype;

      if (first)
      {
        first = false;

        // Create a temporary MPI datatype
        MPI_Datatype tmp_type;
        
        // Define the structure of the pair
        int block_lengths[2] = {1, 1};
        MPI_Aint displacements[2];
        MPI_Datatype block_types[2] = {Mpi_datatype<T1>::value(), Mpi_datatype<T2>::value()};
        
        // Create a pair to calculate memory offsets
        std::pair<T1,T2> dummy_pair;
        
        // Calculate memory displacements using the address of each member
        MPI_Aint base_address;
        MPI_Get_address(&dummy_pair, &base_address);
        MPI_Get_address(&dummy_pair.first, &displacements[0]);
        MPI_Get_address(&dummy_pair.second, &displacements[1]);
        
        // Make relative to the address of the struct itself
        displacements[0] = MPI_Aint_diff(displacements[0], base_address);
        displacements[1] = MPI_Aint_diff(displacements[1], base_address);
        
        // Create the MPI datatype
        MPI_Type_create_struct(2, block_lengths, displacements, block_types, &tmp_type);
        
        // Commit it for use
        MPI_Type_commit(&tmp_type);
        
        // Create a resized type to ensure proper padding/alignment
        MPI_Type_create_resized(tmp_type, 0,  (MPI_Aint)sizeof(dummy_pair), &datatype);
        MPI_Type_commit(&datatype);
        MPI_Type_free(&tmp_type);  // Free the temporary type
      }

      return datatype;
    }
  };
  

  /**
   * @brief Variadic template function to build output string with spaces between arguments. 
   * 
   * All arguments should support << operator.
   * 
   */
  template <typename T, typename... Args>
  void print_log(const T &first, const Args &...args)
  {
    std::ostringstream oss;
    oss << first; // Output the first argument directly
    // Use fold expression to concatenate all arguments with spaces in between
    ((oss << ' ' << args), ...);
    // Output the concatenated string
    std::cout << oss.str() << std::endl;
  }

  /**
   * @brief   Wrapper function on `print_log` that prepends [MPI_rank] to the print output
   * 
   */
  template <typename... Args>
  void print_log_mpi(const int mpi_rank, const Args&... args)
  {
    print_log("[",mpi_rank,"]", args...);
  }

  /**
   * @brief   Builds a string representation from an `std::vector`
   * 
   */
  template <typename T>
  std::string VectorToString(std::vector<T>& vec)
  {
    std::ostringstream output;
    output << "[ ";
    for (auto element : vec)
    {
      output << element << ", ";
    }
    output << "]\n";
    return output.str();
  }

  /**
   * @brief   Builds a string representation from an `std::vector`, only for first `size` number of elements.
   * 
   */
  template <typename T>
  std::string VectorToString(std::vector<T>& vec, size_t size)
  {
    std::ostringstream output;
    output << "[ ";
    for (size_t i = 0; i < size; i++)
    {
      output << vec[i] << ", ";
    }

    output << "]\n";
    return output.str();
  }

  /**
   * @brief Pretty print v1 stats
   * 
   * @param total_time_min 
   * @param total_time_max 
   * @param total_time_avg 
   * @param com_time_min 
   * @param com_time_max 
   * @param com_time_avg 
   * @return std::string 
   */
  std::string FormatStatsV1(uint64_t total_time_min, uint64_t total_time_max, uint64_t total_time_avg,
                          uint64_t com_time_min, uint64_t com_time_max, uint64_t com_time_avg);

  /**
   * @brief Pretty print v2 stats
   * 
   * @param send_size_min 
   * @param send_size_max 
   * @param send_size_avg 
   * @param recv_size_min 
   * @param recv_size_max 
   * @param recv_size_avg 
   * @param sendrecv_size_min 
   * @param sendrecv_size_max 
   * @param sendrecv_size_avg 
   * @return std::string 
   */
  std::string FormatStatsV2(uint64_t send_size_min, uint64_t send_size_max, uint64_t send_size_avg,
                            uint64_t recv_size_min, uint64_t recv_size_max, uint64_t recv_size_avg,
                            uint64_t sendrecv_size_min, uint64_t sendrecv_size_max, uint64_t sendrecv_size_avg);

} // namespace amracut



#endif /* AMRACUT_UTILS_H__ */