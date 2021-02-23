/*+++++++++++++++++ class DoubleParser +++++++++++++++++++ */
#ifndef __doubleparser_h__
#define __doubleparser_h__

#include<math.h>
#include<iostream>

class DoubleParser {
public:
    /* ####  constructors  #### */
    DoubleParser() = default;

    /* ####  destructor  #### */
    ~DoubleParser() = default;

    /* ####  access functions  #### */
    static double valueToCosts(double value, double threshold);

    static double simpleCutOff(double value, double threshold);

private:
    /* ####  member variables  #### */

};

#endif
