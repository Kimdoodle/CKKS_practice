#ifndef FUNCTION_SEAL_H
#define FUNCTION_SEAL_H

#include "SEAL_MM.h"


Ciphertext evaluate_function(vector<double>& poly, Ciphertext& x, int size, string scaleMode, ckks_build& ckks);
Ciphertext sgn_seal(string mode, string scaleMode, Ciphertext& x, vector<double> poly, vector<double>& input, int d, int pre, ckks_build& ckks);
Ciphertext abs_seal(string mode, string scaleMode, Ciphertext& x, vector<double>& poly, vector<double>& input, int d, int pre, ckks_build& ckks);
Ciphertext max_seal(string mode, string scaleMode, Ciphertext& x, vector<double>& poly, vector<double>& input, int d, int pre, ckks_build& ckks);

#endif // FUNCTION_SEAL_H