#pragma once

#include <string>

#include "cluster_editing/definitions.h"

namespace cluster_editing {
namespace io {

Graph readGraphFile(const std::string& filename);

}
} // namespace cluster_editing::io