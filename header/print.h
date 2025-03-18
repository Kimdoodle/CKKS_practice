#ifndef PRINT_H
#define PRINT_H

#include "SEAL_VS.h"

void printVector(vector<double>& coeffs, bool asFunction, int pre = 4);
void print_parameters(const SEALContext& context);
void printStep(Ciphertext& y, vector<double>& real_result, vector<double>& poly, int iter, int size, ckks_build& ckks, int pre);
void printResult(Ciphertext& y, vector<double>& real_result, int size, ckks_build& ckks, int pre);

#endif // PRINT_H
