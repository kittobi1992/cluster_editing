/*+++++++++++++++++ class ProblemInstanceException +++++++++++++++++++ */

#include <Exceptions/probleminstanceexception.h>

ProblemInstanceException::ProblemInstanceException(std::string message) : _message(std::move(message)) {}

std::string ProblemInstanceException::getMessage() { 
	return _message; 
}
