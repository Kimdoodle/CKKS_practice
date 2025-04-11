#ifndef POLYMATH_H
#define POLYMATH_H

#include "SEAL_MM.h"

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

vector<vector<double>> make_matrix(int inputsize, const string& type);
vector<vector<double>> type_matrix(int originSize, const string& type);
vector<vector<double>> pad_matrix(vector<vector<double>> U, int d);
vector<vector<double>> ipad_matrix(vector<vector<double>> U, int d);
vector<double> flatten_matrix(vector<vector<double>> U);
vector<vector<double>> unflatten_matrix(vector<double> U, int colsize);
vector<vector<double>> diagonal_matrix(vector<vector<double>> U);


#endif // POLYMATH_H
