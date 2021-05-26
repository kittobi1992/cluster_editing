#pragma once

#include <string>

#include "cluster_editing/definitions.h"

namespace cluster_editing::io {

Graph readGraphFile();
Graph readGraphFile(const std::string& filename);

} // namespace cluster_editing::io