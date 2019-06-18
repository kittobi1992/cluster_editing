#include <doubleparser.h>

double DoubleParser::valueToCosts(double value, double threshold)
{
	return value-threshold;
}

double DoubleParser::simpleCutOff(double value, double threshold)
{
	return (value < threshold) ? -1.0 : 1.0;
}
