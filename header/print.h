#ifndef PRINT_H
#define PRINT_H

#include "SEAL_MM.h"

void printVector(vector<double>& coeffs, bool asFunction, int pre = 4);
void printMatrix(vector<vector<double>> U, int rowsize=0);
void print_parameters(const SEALContext& context);
void printStep(vector<double>& realValue, vector<double>& poly, vector<double>& fnDec, vector<double> raw_inputs, Ciphertext& y, ckks_build& ckks, string function, int i);

#endif // PRINT_H
