#include <algorithm>
#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>
#include <filesystem>
#include <iostream>

#include "cluster_editing/io/draw_graph.h"

namespace po = boost::program_options;

int main(int argc, char* argv[]) {
  // command line options
  po::options_description options("Required Options");
  std::filesystem::path graph_filename;
  std::filesystem::path coords_filename;
  std::filesystem::path partition_filename;
  std::filesystem::path output_filename;
  options.add_options()
      ("help", "show help message")
      ("graph,g",
       po::value<std::filesystem::path>(&graph_filename)
       ->value_name("<path>")
       ->required(),
       "graph filename")
      ("coordinates,c",
       po::value<std::filesystem::path>(&coords_filename)
       ->value_name("<path>")
       ->required(),
       "coordinates filename")
      ("partition,p",
       po::value<std::filesystem::path>(&partition_filename)
       ->value_name("<path>")
       ->required(),
       "partition filename")
      ("output,o",
       po::value<std::filesystem::path>(&output_filename)
       ->value_name("<path>")
       ->required(),
       "output filename");

  po::variables_map cmd_vm;
  po::store(po::parse_command_line(argc, argv, options), cmd_vm);

  if (cmd_vm.count("help") != 0 || argc == 1) {
    std::cout << options << std::endl;
    exit(0);
  }
  
  po::notify(cmd_vm);

  // drawing the graph
  draw_graph(graph_filename, coords_filename, partition_filename, output_filename);

  return 0;
}
