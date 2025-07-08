#include <sstream>
#include <algorithm>
#include <iomanip>

#include "utils.hpp"

std::string amracut::FormatStatsV1(uint64_t total_time_min, 
                                   uint64_t total_time_max, 
                                   uint64_t total_time_avg)
{
  size_t min_col_width = std::max({std::string("min").length(),
                                   std::to_string(total_time_min).length()});

  size_t max_col_width = std::max({std::string("max").length(),
                                   std::to_string(total_time_max).length()});

  size_t avg_col_width = std::max({std::string("avg").length(),
                                   std::to_string(total_time_avg).length()});

  size_t label_col_width = std::string("total time (us)").length();

  const size_t padding = 2;
  min_col_width += padding;
  max_col_width += padding;
  avg_col_width += padding;

  std::ostringstream output;

  // Print header
  output << std::setw(label_col_width) << " " // Space for row labels
         << std::setw(min_col_width) << "min"
         << std::setw(max_col_width) << "max"
         << std::setw(avg_col_width) << "avg" << '\n';

  // Print total time values
  output << std::setw(label_col_width) << "total time (us)"
         << std::setw(min_col_width) << total_time_min
         << std::setw(max_col_width) << total_time_max
         << std::setw(avg_col_width) << total_time_avg << '\n';


  return output.str();
}

