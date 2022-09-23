#include "filters.h"

void medianFilter(std::vector<double> &outData, std::vector<double> &inData, const int &windowSize)
{
	/*! Function describing the median filter
	* Inputs:
	*		inData		- processed data, type std::vector<double>
	*		windowSize		- size of filtering window, type int
	* Outputs:
	*		outData		- filtered data, type type std::vector<double>
	*/

	assert(!inData.empty()); //Check processed data
	std::vector<double> window(windowSize);
	std::vector<double>::iterator iterator;
	outData.clear(); // clear output data vector
	for (iterator = inData.begin(); iterator != inData.begin() + windowSize / 2; iterator++) outData.push_back(*iterator);
	for (iterator = inData.begin() + windowSize / 2; iterator != inData.end() - windowSize / 2; iterator++)
	{
		window.assign(iterator - windowSize / 2, iterator + windowSize / 2 + 1);
		std::nth_element(window.begin(), window.begin() + windowSize / 2, window.end());
		outData.push_back(*(window.begin() + windowSize / 2));
	}
	for (iterator = inData.end() - windowSize / 2; iterator != inData.end(); iterator++) outData.push_back(*iterator);
}

void butterworthFilter(std::vector<double> &outData, const std::vector<double> &inData)
{
	/*! Function describing the Butterworth filter 4 order (one section)
	* Inputs:
	*		inData		- input signal
	* Outputs:
	*		outData		- the Butterworth filtered signal
	*/

	assert(!inData.empty()); //Check parameters
	outData.clear(); // clear output vector
	double NUM[] = {6.30129962118151e-06, 1.89038988635445e-05, 1.89038988635445e-05, 6.30129962118151e-06};
    double DEN[] = {1, -2.92520452880893, 2.85318010843919, -0.927925169233285};

    double x, y0, y1 = 0, y2 = 0, y3 = 0;
    double b0= NUM[0]; double b1= NUM[1]; double b2= NUM[2]; double b3= NUM[3]; 
	double a0= DEN[0]; double a1= DEN[1]; double a2= DEN[2]; double a3= DEN[3];

	for (size_t i = 0; i < inData.size(); i++) // for each element of a input vector do
    {
        x = inData.at(i); 
		// coefficients calculation
        y0 = x * b0 + y1;
        y1 = b1 * x - a1 * y0 + y2;
        y2 = b2 * x - a2 * y0 + y3;
        y3 = b3 * x - a3 * y0;
		// calculations done
        outData.push_back(y0 * a0); // computing filtered value and add in outData
    }
}

void savitzkyGolayFilter(std::vector<double> &outData, std::vector<double> &inData, const int &windowSize, const int &polyDegree)
{
	/*! Function describing the Savitzky-Golay filter
	* Inputs:
	*		inData		- input signal
	*		windowSize		- size of window or polynomial order
	* Outputs:
	*		outData		- filtered signal
	*/
	assert (!inData.empty ()); //Check parameters
	assert ((windowSize == 5 && polyDegree == 2) || (windowSize == 7 && polyDegree == 2) || (windowSize == 7 && polyDegree == 4));
	double *polyCoeff = new double[windowSize];
	double normCoeff;
	if ( windowSize == 5 && polyDegree == 2 )
	{
		polyCoeff[0] = -3.0;
		polyCoeff[1] = 12.0;
		polyCoeff[2] = 17.0;
		polyCoeff[3] = 12.0;
		polyCoeff[4] = -3.0;
		normCoeff = 35.0;
	}
	else
	if ( windowSize == 7 && polyDegree == 2 )
	{
		polyCoeff[0] = -2.0;
		polyCoeff[1] = 3.0;
		polyCoeff[2] = 6.0;
		polyCoeff[3] = 7.0;
		polyCoeff[4] = 6.0;
		polyCoeff[5] = 3.0;
		polyCoeff[6] = -2.0;
		normCoeff = 21.0;
	}
	else
	if ( windowSize == 7 && polyDegree == 4 )
	{
		polyCoeff[0] = 5.0;
		polyCoeff[1] = -30.0;
		polyCoeff[2] = 75.0;
		polyCoeff[3] = 131.0;
		polyCoeff[4] = 75.0;
		polyCoeff[5] = -30.0;
		polyCoeff[6] = 5.0;
		normCoeff = 231.0;
	}


	outData.clear (); // clear output vector
	std::vector<double> tmpArray (windowSize);
	std::vector<double>::iterator iter, tmpIter;
	// доопределяем циклически краевые значения
	std::vector<double> tmpSource;
	tmpSource.resize (inData.size () + windowSize - 1);
	int k = 0;
	for ( size_t i = 0; i < inData.size () + windowSize - 1; i++ )
	{
		if ( i < (size_t) windowSize / 2 || i > inData.size () + windowSize / 2 - 1 )
			tmpSource.at (i) = 1.0;
		else
		{
			tmpSource.at (i) = inData.at (k);
			k++;
		}
	}
	// сглаживаем скользящим окном
	for ( iter = tmpSource.begin () + windowSize / 2; iter != tmpSource.end () - windowSize / 2; iter++ )
	{
		tmpArray.assign (iter - windowSize / 2, iter + windowSize / 2 + 1);
		double sum = 0.0;
		int k = 0;
		for ( tmpIter = tmpArray.begin (); tmpIter != tmpArray.end (); tmpIter++ )
		{
			sum = sum + *tmpIter * polyCoeff[k];
			k++;
		}
		sum = sum / normCoeff;
		outData.push_back (sum);
	}
	delete[] polyCoeff;
}

void grubbsFilter (std::vector<double> &outData, std::vector<double> &inData, const double &alpha)
{
	double threshold;
	if ( alpha == 0.1 )
		threshold = 3, 53;
	else if ( alpha == 0.05 )
		threshold = 3.70;
	else if ( alpha == 0.01 )
		threshold = 4.11;
	double sum = 0, sum2 = 0;
	std::vector<double>::iterator iter;
	for ( iter = inData.begin(); iter != inData.end (); iter++ )
	{
		sum += *iter;
		sum2 += pow ((*iter), 2);
	}
	double mean = sum / inData.size ();
	double mean2 = sum2 / inData.size ();
	double sigma = mean2 - pow (mean, 2);
	for ( iter = inData.begin (); iter < inData.end (); iter += 1 ) // for each element of input vector do
	{
		if ( ((*iter) - mean) / sigma <= threshold )
			outData.push_back (*iter);
	}
}