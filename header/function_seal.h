#pragma once

#include "SEAL_VS.h"


Ciphertext evaluate_function(vector<double>& poly, Ciphertext& x, int size, string scaleMode, ckks_build& ckks);
Ciphertext sgn_seal(string mode, string scaleMode, Ciphertext& x, vector<double> poly, vector<double>& input, int d, int pre, ckks_build& ckks);
Ciphertext abs_seal(string mode, string scaleMode, Ciphertext& x, vector<double>& poly, vector<double>& input, int d, int pre, ckks_build& ckks);
Ciphertext max_seal(string mode, string scaleMode, Ciphertext& x, vector<double>& poly, vector<double>& input, int d, int pre, ckks_build& ckks);

void newton_seal(double x, vector<double> x0, int iter, string printmode, ckks_build& ckks);
void goldschmidt_seal(double x, vector<double> x0, int iter, string printmode, ckks_build& ckks);