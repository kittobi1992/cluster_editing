/*+++++++++++++++++ class GraphException +++++++++++++++++++ */
#ifndef __graphException_h__
#define __graphException_h__

#include <string>

/*
Simply implements a container for a string message,
which I call Exception in generell.
This idea copies somehow the concept of Exceptions
from Java.
*/

class GraphException : std::exception {
public:
	// Constructor
	explicit GraphException(std::string msg);

	// access function
	std::string getMessage();

private:
	// error message
	std::string message;
};

#endif