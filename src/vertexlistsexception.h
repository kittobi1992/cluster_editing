/*+++++++++++++++++ class VertexListsException +++++++++++++++++++ */
#ifndef __vertexlistsexception_h__
#define __vertexlistsexception_h__

#include <string>

/*
Simply implements a container for a string message,
which I call Exception in generell.
This idea copies somehow the concept of Exceptions
from Java.
*/

class VertexListsException {
public:
	// Constructor
	VertexListsException(std::string message);

	// Destructor
	~VertexListsException() {};

	// access function
	std::string getMessage();

private:
	// error message
	std::string _message;
};

#endif
