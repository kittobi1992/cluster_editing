/*+++++++++++++++++ class GraphException +++++++++++++++++++ */

#include <graphexception.h>

// Constructor
GraphException::GraphException(std::string msg) : message(msg) { }

// Destructor
GraphException::~GraphException() { }

// return error message
std::string GraphException::getMessage()
{
	return this->message;
};
