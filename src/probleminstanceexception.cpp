/*+++++++++++++++++ class ProblemInstanceException +++++++++++++++++++ */

#include <probleminstanceexception.h>

ProblemInstanceException::ProblemInstanceException(std::string message) : _message(message) {
};

std::string ProblemInstanceException::getMessage() { 
	return _message; 
};
