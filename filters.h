#include <vector>
#include <algorithm>
#include <assert.h>

void medianFilter(std::vector<double> &outData, std::vector<double> &inData, const int &windowSize);
void butterworthFilter (std::vector<double> &outputData, const std::vector<double> &inputData);
void savitzkyGolayFilter (std::vector<double> &outData, std::vector<double> &inData, const int &windowSize, const int &polynomDegree);
void grubbsFilter (std::vector<double> &outData, std::vector<double> &inData, const double &alpha);