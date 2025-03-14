#ifndef FUNCTION_SEAL_H
#define FUNCTION_SEAL_H

#include "SEAL_VS.h"

Ciphertext evaluate_function(vector<double>& poly, Ciphertext& x, ckks_build& ckks);
Ciphertext sgn(string mode, Ciphertext& x, vector<double>& poly, int d, ckks_build& ckks);
Ciphertext abs_seal(Ciphertext& input, vector<double>& poly, int d, ckks_build& ckks);
Ciphertext max_seal(vector<double>& input, vector<double>& poly, int d, ckks_build& ckks);

#endif // FUNCTION_SEAL_H