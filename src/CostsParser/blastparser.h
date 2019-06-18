/*+++++++++++++++++ class BlastParser +++++++++++++++++++ */
#ifndef __blastparser_h__
#define __blastparser_h__

#include<math.h>
#include<iostream>

class BlastParser {
public:
	/* ####  constructors  #### */
	BlastParser();
	/* ####  destructor  #### */
	~BlastParser();
	
	/* ####  access functions  #### */
	static double exponentToCosts(double blast_value, double threshold);
	static double simpleCutOff(double blast_value, double threshold);

private:
	/* ####  member variables  #### */
	
};

#endif
