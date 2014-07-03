#include "AquaTopBGGenerator.h"
#include <math.h>

AquaTopBGGenerator::AquaTopBGGenerator()
{
}

AquaTopBGGenerator::~AquaTopBGGenerator()
{
	
}

// Fitting target: lowest sum of squared absolute error
// Fitting target value = 0.747158511837
//z = a + by + cx + dxy + f(x^2) + g(x^2)y
double AquaTopBGGenerator::GetZValue(double x_in, double y_in)
{
	double temp;
	temp = 0.0;

	// coefficients
	double a = 1.8564803894488152E+03;
	double b = 5.2563093116817700E-02;
	double c = 1.2853574330037623E-01;
	double d = -2.7047119933664614E-04;
	double f = -1.9683753259302163E-04;
	double g = 4.7046643987525005E-07;

	temp += a;
	temp += b * y_in;
	temp += c * x_in;
	temp += d * x_in * y_in;
	temp += f * pow(x_in, 2.0);
	temp += g * pow(x_in, 2.0) * y_in;
	return temp;
}

//z = a + bx + cy
double AquaTopBGGenerator::GetZValueLin(double x_in, double y_in)
{
	double temp;
	temp = 0.0;

	// coefficients
	double a = 1.8706125660880882E+03;
	double b = 1.9435859943767713E-02;
	double c = 1.6705938888233653E-02;

	temp = a + b * x_in + c * y_in;
	return temp;
}