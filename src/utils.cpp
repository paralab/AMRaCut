#include <sstream>
#include <algorithm>
#include <iomanip>

#include "utils.hpp"

std::string amracut::FormatStatsV1(uint64_t total_time_min, uint64_t total_time_max, uint64_t total_time_avg,
                                    uint64_t com_time_min, uint64_t com_time_max, uint64_t com_time_avg)
{
  size_t min_col_width = std::max({std::string("min").length(),
                                   std::to_string(total_time_min).length(),
                                   std::to_string(com_time_min).length()});

  size_t max_col_width = std::max({std::string("max").length(),
                                   std::to_string(total_time_max).length(),
                                   std::to_string(com_time_max).length()});

  size_t avg_col_width = std::max({std::string("avg").length(),
                                   std::to_string(total_time_avg).length(),
                                   std::to_string(com_time_avg).length()});

  size_t label_col_width = std::max(std::string("total time (us)").length(),
                                    std::string("com time (us)").length());

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

  // Print com time values
  output << std::setw(label_col_width) << "com time (us)"
         << std::setw(min_col_width) << com_time_min
         << std::setw(max_col_width) << com_time_max
         << std::setw(avg_col_width) << com_time_avg;

  return output.str();
}

std::string amracut::FormatStatsV2(uint64_t send_size_min, uint64_t send_size_max, uint64_t send_size_avg,
                                    uint64_t recv_size_min, uint64_t recv_size_max, uint64_t recv_size_avg,
                                    uint64_t sendrecv_size_min, uint64_t sendrecv_size_max, uint64_t sendrecv_size_avg)
{
  size_t min_col_width = std::max({std::string("min").length(),
                                   std::to_string(send_size_min).length(),
                                   std::to_string(recv_size_min).length(),
                                   std::to_string(sendrecv_size_min).length()});

  size_t max_col_width = std::max({std::string("max").length(),
                                   std::to_string(send_size_max).length(),
                                   std::to_string(recv_size_max).length(),
                                   std::to_string(sendrecv_size_max).length()});

  size_t avg_col_width = std::max({std::string("avg").length(),
                                   std::to_string(send_size_avg).length(),
                                   std::to_string(recv_size_avg).length(),
                                   std::to_string(sendrecv_size_avg).length()});

  size_t label_col_width = std::max({std::string("sent (bytes)").length(),
                                    std::string("received (bytes)").length(),
                                    std::string("total (bytes)").length()});

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

  // Print total send size values
  output << std::setw(label_col_width) << "sent (bytes)"
         << std::setw(min_col_width) << send_size_min
         << std::setw(max_col_width) << send_size_max
         << std::setw(avg_col_width) << send_size_avg << '\n';

  // Print recev size values
  output << std::setw(label_col_width) << "received (bytes)"
         << std::setw(min_col_width) << recv_size_min
         << std::setw(max_col_width) << recv_size_max
         << std::setw(avg_col_width) << recv_size_avg << '\n';

  // Print total send+recv size values
  output << std::setw(label_col_width) << "total (bytes)"
         << std::setw(min_col_width) << sendrecv_size_min
         << std::setw(max_col_width) << sendrecv_size_max
         << std::setw(avg_col_width) << sendrecv_size_avg;

  return output.str();
}