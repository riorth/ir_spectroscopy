#define _CRT_SECURE_NO_WARNINGS // отключение предупреждений
#include <iostream> // включаем функции для организации ввода и вывода
#include <sstream> // включаем функции для работы со строками
#include <fstream> // влючаем функции работы с файлами
#include <vector> // библиотека для работы с векторами
#include <algorithm> // библиотека для работы с дополнительными функциями
#include <assert.h> // для проверки входных данных

typedef std::vector<double> vector;
typedef std::vector<std::vector<double>> matrix;

std::vector<double>  medianFilter (std::vector<double> &spectrum_in, const int &windowSize)
{
	// Медианный фильтр. На вход подается вектор со значениями спектра и параметр «размер окна», на выходе получается вектор с отфильтрованными значениями
	std::vector<double> result;
	assert (!spectrum_in.empty ()); // Проверяем чтобы входной вектор был не пуст
	result.clear (); // отчищаем вектор в который записывается результат
	std::vector<double> tmpArray (windowSize);
	std::vector<double>::iterator iterator;
	for ( iterator = spectrum_in.begin (); iterator != spectrum_in.begin () + windowSize / 2; iterator++ ) result.push_back (*iterator);
	for ( iterator = spectrum_in.begin () + windowSize / 2; iterator != spectrum_in.end () - windowSize / 2; iterator++ )
	{
		tmpArray.assign (iterator - windowSize / 2, iterator + windowSize / 2 + 1);
		std::nth_element (tmpArray.begin (), tmpArray.begin () + windowSize / 2, tmpArray.end ());
		result.push_back (*(tmpArray.begin () + windowSize / 2));
	}
	for ( iterator = spectrum_in.end () - windowSize / 2; iterator != spectrum_in.end (); iterator++ ) result.push_back (*iterator);
	return result;
}

std::vector<double> butterworthFilter (const std::vector<double> &spectrum_in)
{
	// Фильтр Баттерворта второго порядка. На вход подается вектор со значениями спектра.
	assert (!spectrum_in.empty ()); // Проверяем чтобы входной вектор был не пуст
	std::vector<double> result;
	result.clear (); // отчищаем вектор в который записывается результат
	double GAIN = 1.992200736e+05;
	double x, y0, x0, x1 = 0, x2 = 0, y1 = 0, y2 = 0, y3 = 0;
	double a0 = -0.9936731239; double a1 = 1.9936530456;
	for ( size_t i = 0; i < spectrum_in.size (); i++ ) // for each element of a input vector do
	{
		x = spectrum_in.at (i);
		x0 = x1;
		x1 = x2;
		x2 = x / GAIN;
		y0 = y1;
		y1 = y2;
		y2 = x0 + x2 + 2 * x1 + a0*y1 + a1*y1;
		result.push_back (y2); // computing filtered value and add in outputData
	}
	return result;
}

std::vector<double> savitzkyGolayFilter (std::vector<double> &spectrum_in, const int &windowSize, const int &polyDegree)
{
	// Фильтр Cавицкого Голея. На вход подается вектор со значениями спектра, размер окна и степень полинома
	assert (!spectrum_in.empty ()); // Проверяем чтобы входной вектор был не пуст
	assert ((windowSize == 5 && polyDegree == 2) || (windowSize == 7 && polyDegree == 2) || (windowSize == 7 && polyDegree == 4));
	std::vector<double> result;
	result.clear (); // отчищаем вектор в который записывается результат
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

	std::vector<double> tmpArray (windowSize);
	std::vector<double>::iterator iterator, tmpIterator;
	// доопределяем циклически краевые значения
	std::vector<double> tmpSource;
	tmpSource.resize (spectrum_in.size () + windowSize - 1);
	int k = 0;
	for ( size_t i = 0; i < spectrum_in.size () + windowSize - 1; i++ )
	{
		if ( i < (size_t) windowSize / 2 || i > spectrum_in.size () + windowSize / 2 - 1 )
			tmpSource.at (i) = 1.0;
		else
		{
			tmpSource.at (i) = spectrum_in.at (k);
			k++;
		}
	}
	// сглаживаем скользящим окном
	for ( iterator = tmpSource.begin () + windowSize / 2; iterator != tmpSource.end () - windowSize / 2; iterator++ )
	{
		tmpArray.assign (iterator - windowSize / 2, iterator + windowSize / 2 + 1);
		double sum = 0.0;
		int k = 0;
		for ( tmpIterator = tmpArray.begin (); tmpIterator != tmpArray.end (); tmpIterator++ )
		{
			sum = sum + *tmpIterator * polyCoeff[k];
			k++;
		}
		sum = sum / normCoeff;
		result.push_back (sum);
	}
	delete[] polyCoeff;
	return result;
}

std::vector<double> grubbsFilter (std::vector<double> &spectrum_in, const double &alpha)
{
	// Критерий Граббса. На вход подается вектор со значениями спектра и задаваемы уровень значимости.
	assert (!spectrum_in.empty ()); // Проверяем чтобы входной вектор был не пуст
	std::vector<double> result;
	result.clear (); // отчищаем вектор в который записывается результат
	double threshold;
	if ( alpha == 0.1 )
		threshold = 3.53;
	else if ( alpha == 0.05 )
		threshold = 3.70;
	else if ( alpha == 0.01 )
		threshold = 4.11;
	double sum = 0, sum2 = 0;
	std::vector<double>::iterator iter;
	for ( iter = spectrum_in.begin (); iter != spectrum_in.end (); iter++ )
	{
		sum += *iter;
		sum2 += pow ((*iter), 2);
	}
	double mean = sum / spectrum_in.size ();
	double mean2 = sum2 / spectrum_in.size ();
	double sigma = mean2 - pow (mean, 2);
	for ( iter = spectrum_in.begin (); iter < spectrum_in.end (); iter += 1 ) // for each element of input vector do
	{
		if ( ((*iter) - mean) / sigma <= threshold )
			result.push_back (*iter);
	}
	return result;
}

int write_file (std::string outPath, std::vector<double> &spectrum)
{
	// фукнция осуществляет запись в спектра в файл
	// входные параметры – путь к файлу в который будет записан спектр и спектр в формате std::vector<double>
	FILE *file;
	file = fopen (outPath.c_str (), "a");
	std::ofstream ofstream_out (outPath);
	for ( size_t i = 0; i < spectrum.size (); i++ )
		ofstream_out << spectrum[i] << ", ";
	fclose (file);
	return 0;
}

void read_file (std::string inPath, std::vector<double> &spectrum)
{
	// фукнция осуществляет чтение спектра из файла
	// входные параметры – путь к файлу из которого будет прочитан спектр и вектор std::vector<double> в который будет записан спектр

	std::ifstream  ifstream_file (inPath);
	std::string line;
	std::getline (ifstream_file, line);
	std::stringstream stringstream_line (line);
	std::string  cell;

	while ( std::getline (stringstream_line, cell, ',') )
	{
		spectrum.push_back (atof (cell.c_str ()));
	};
};
bool exists_file (std::string file_name)
{
	// функция проверяет существования файла, входным параметром является путь к файлу
	return std::ifstream (file_name).good ();
};
// ниже две функции необходимые для конвертирования строки в числа и обратно
template <typename Type, typename String>
Type StrToType (const String &string)
{
	Type variable;
	std::stringstream stringStream (string);
	stringStream >> variable;
	return variable;
}

template <typename Type>
std::string TypeToStr (const Type &variable)
{
	std::stringstream stringStream;
	stringStream << variable;
	return stringStream.str ();
}
// основная функция программы
int main (int argc, char * argv[])
{
	std::string inPath = "";
	std::string outPath = "";
	std::string type_filter = "";
	int window_size (0);
	int poly_degree (0);
	double alpha (0);
	if ( argc < 3 )
	{
		std::cout << "\n   Help to start:\n";
		std::cout << "   There is no input parameters\n";
		std::cout << "   You can take advantage of \? or -help, to call help\n";
		std::cout << "   -in  <path to file>  - the file path with the spectrum processing\n";
		return 1;
	}
	else
	{
		int i = 1;
		while (i < argc)
		{
			if ( argv[i] == std::string ("-help") || argv[i] == std::string ("/?") )
			{
				std::cout << "\n   Help to start:\n";
				std::cout << "   -in  <path to file>  - the file path with the spectrum processing\n";
				std::cout << "   -out <path to file>  - the file path with the spectrum after processing\n";
				std::cout << "   -type  <type of filter>  - medianFilter, butterworthFilter, savitzkyGolayFilter, grubbsFilter\n";
				std::cout << "   -params <param1, param2>		 - medianFilter(windowSize), savitzkyGolayFilter(windowSize, polyDegree), butterworthFilter(), grubbsFilter(alpha)\n";
				return 1;
			}
				if (argv[i] == std::string("-in"))
				{
					inPath = argv[i + 1];
					if ( !exists_file(inPath) )
					{
						std::cout << "The path with the spectrum is invalid\n";
						return 1;
					}
					i += 2;

				}
				else if (argv[i] == std::string("-out"))
				{
					outPath = argv[i + 1];
					if ( outPath != "" )
					{
						std::ifstream fin (outPath, std::ios_base::out | std::ios_base::trunc);
						if ( !fin.is_open () )
						{
							std::cout << "The file path with the spectrum after processing  is invalid\n";
							return 1;
						}
						else
							fin.close ();
					}
					i += 2;
				}
				else if ( argv[i] == std::string ("-type") )
				{
					type_filter = argv[i + 1];
					if ( type_filter != "medianFilter" && type_filter != "savitzkyGolayFilter" && type_filter != "butterworthFilter" && type_filter != "grubbsFilter" )
					{
						std::cout << "The type of filter is invalid\n";
						return 1;
					}
					i += 2;
				}
				else if ( argv[i] == std::string ("-params") )
				{
					if ( type_filter == "medianFilter" )
					{
						if ( argc <= i + 1 )
						{
							std::cout << "The params of filter is invalid\n";
							return 1;
						}							
						else
							window_size = StrToType<int> (argv[i + 1]);						
					}
					else if ( type_filter == "savitzkyGolayFilter" )
					{
						if ( argc <= i + 1 )
						{
							std::cout << "The params of filter is invalid\n";
							return 1;
						}
						else
						{
							std::string params = argv[i + 1];
							window_size = StrToType<int> (params.substr (0, params.find_first_of (",")));
							poly_degree = StrToType<int> (params.substr (params.find_first_of (","), params.size ()));
						}
					}
					else if ( type_filter == "grubbsFilter" )
					{
						if ( argc <= i + 1 )
						{
							std::cout << "The params of filter is invalid\n";
							return 1;
						}
						else
						{
							alpha = StrToType<double> (argv[i + 1]);
						}
					}
					else
					{
						std::cout << "No valid argument\n";
						return 1;
					}
					i += 2;
				}
				else
				{
					std::cout << "No valid argument\n";
					return 1;
				}
		}

	}
	std::vector<double> spectrum_in;
	std::vector<double> spectrum_out;
	if ( inPath != "" && outPath != "" && type_filter != "" )
	{
		read_file (inPath, spectrum_in);
		if ( type_filter == "medianFilter" && window_size > 0 && window_size < static_cast<int>(spectrum_in.size ()) )
		{
			spectrum_out = medianFilter (spectrum_in, window_size);
		}
		else if ( type_filter == "butterworthFilter" )
		{
			spectrum_out = butterworthFilter (spectrum_in);
		}
		else if ( type_filter == "savitzkyGolayFilter" && ((window_size == 5 && poly_degree == 2) || (window_size == 7 && poly_degree == 2) || (window_size == 7 && poly_degree == 4)))
		{
			spectrum_out = savitzkyGolayFilter (spectrum_in, window_size, poly_degree);
		}
		else if ( type_filter == "grubbsFilter" && ((alpha == 0.1 || alpha == 0.01 || alpha == 0.05)))
		{
			spectrum_out = grubbsFilter (spectrum_in, alpha);
		}
		else
		{
			std::cout << "incorrect parameters or input data\n";
			return 1;
		}
	};
	write_file (outPath, spectrum_out);
	std::cout << "Filtering completed successfully" << std::endl;
}