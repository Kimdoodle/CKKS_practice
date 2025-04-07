#ifndef POLYMATH_H
#define POLYMATH_H

#include "SEAL_VS.h"

int factorial(int a, int b=0);
double log2(double x, double base);
vector<double> differentiate(vector<double> poly);
vector<double> sample_data(double min, double max, double epsilon = 0.0, int iter = 1);
vector<double> duplicate_vector(double input, int size);
vector<vector<double>> createToeplitzMatrix(const vector<double>& coeffs, int result_size);
vector<double> multiplyMatrixVector(const vector<vector<double>>& matrix, const vector<double>& vec);
vector<double> multPlainPolynomial(vector<double>& v, double scalar);
vector<double> multPolynomial(const vector<double>& a, const vector<double>& b);
vector<double> powerPolynomial(const vector<double>& poly, int exponent);
double polyEvaluate(const vector<double>& poly, double input);
double polypolyEvaluate(const vector<double>& poly, double input, int d);
vector<double> polypolyEvaluate(const vector<double>& poly, vector<double>& input);
vector<double> calculatePoly(const vector<double>& x, const vector<double>& y);

#endif // POLYMATH_H
