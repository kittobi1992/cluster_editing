#pragma once

#include <string>

#include "cluster_editing/definitions.h"

namespace cluster_editing::io {

Graph readGraphFile(const std::string& filename);

Graph readGraphFromStandardInput();

} // namespace cluster_editing::io