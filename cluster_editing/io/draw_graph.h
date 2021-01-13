#pragma once

#include "cluster_editing/definitions.h"
#include <filesystem>
#include <vector>

using namespace cluster_editing;

/// Draw a given graph with given coordinates and partition.
/// \param graph_file_in Graph file in the PACE-format.
/// \param coords_file_in File containing vertex coordinates. Each line should
///          consist of three numbers; first comes the vertex id (starting with
///          1 as in the graph file) and then the x and y coordinate.
/// \param partition_file_in File containing the vertex partition. Each line
///          should consist of two integers; first comes the vertex id (starting
///          with 1 as in the graph file) and then the partition id.
/// \param ipe_file_out File that should store the resulting ipe file.
void draw_graph(const std::filesystem::path& graph_file_in,
                const std::filesystem::path& coords_file_in,
                const std::filesystem::path& partition_file_in,
                const std::filesystem::path& ipe_file_out);

/// As the other draw_graph but with the internal data structures instead of
/// reading everything from file.
void draw_graph(const Graph& G,
                const std::vector<std::vector<double>>& coords,
                const std::vector<unsigned>& partition,
                const std::filesystem::path& ipe_file_out);

std::vector<std::vector<double>> read_coords(
    const std::filesystem::path& coords_file_in, unsigned num_nodes);

std::vector<unsigned> read_partition(
    const std::filesystem::path& partition_file_in, unsigned num_nodes);
