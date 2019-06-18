/*+++++++++++++++++ class VertexListsException +++++++++++++++++++ */

#include <vertexlistsexception.h>

/* #################  Constructor  ##################### */

VertexListsException::VertexListsException(std::string message) : _message(message) {
}

// access function
std::string VertexListsException::getMessage() { 
	return _message;
}
