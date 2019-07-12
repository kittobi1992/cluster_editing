/*+++++++++++++++++ class GraphException +++++++++++++++++++ */

#include <Exceptions/graphexception.h>

// Constructor
GraphException::GraphException(std::string msg) : message(std::move(msg)) {}

// return error message
std::string GraphException::getMessage() {
    return message;
}
