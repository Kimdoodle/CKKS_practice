#pragma once

#include "SEAL_VS.h"

void printVector(vector<double>& coeffs, bool asFunction, int pre = 4);
void print_parameters(const SEALContext& context);
void printStep(vector<double>& realValue, vector<double>& poly, vector<double>& fnDec, vector<double> raw_inputs, Ciphertext& y, ckks_build& ckks, string function, int i);

void debug_print(string message, string mode);
string toScientific(double num, int precision);
void printVector_10eform(const vector<double>& coeffs, bool asFunction, int pre=4);
