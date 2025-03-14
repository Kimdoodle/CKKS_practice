#ifndef POLYMATH_H
#define POLYMATH_H

#include "SEAL_VS.h"

int factorial(int a);
double log2(double x, double base);
vector<vector<double>> createToeplitzMatrix(const vector<double>& coeffs, int result_size);
vector<double> multiplyMatrixVector(const vector<vector<double>>& matrix, const vector<double>& vec);
vector<double> multPolynomial(const vector<double>& a, const vector<double>& b);
vector<double> powerPolynomial(const vector<double>& poly, int exponent);
double polyEvaluate(const vector<double>& poly, double input);
double polypolyEvaluate(const vector<double>& poly, double input, int d);

#endif // POLYMATH_H
