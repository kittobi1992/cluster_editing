#include <blastparser.h>

double BlastParser::exponentToCosts(double blast_value, double threshold) {
    double log1 = log10(blast_value + 1e-300);
    double result = (log10(threshold) - log1);
    return (result == 0.0) ? 1e-300 : result;
    //return result;
}

double BlastParser::simpleCutOff(double blast_value, double threshold) {
    return (blast_value < threshold) ? 1.0 : -1.0;
}
