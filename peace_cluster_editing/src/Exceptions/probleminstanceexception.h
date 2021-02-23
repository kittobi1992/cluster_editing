/*+++++++++++++++++ class ProblemInstanceException +++++++++++++++++++ */
#ifndef __probleminstanceexception_h__
#define __probleminstanceexception_h__

#include <string>

/*
Simply implements a container for a string message,
which I call Exception in generell.
This idea copies somehow the concept of Exceptions
from Java.
*/

class ProblemInstanceException : std::exception {
public:
    /* ####  Constructor  #### */
    explicit ProblemInstanceException(std::string message);

    /* ####  access function  #### */
    std::string getMessage();

private:
    std::string _message;
};


#endif
