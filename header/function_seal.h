#ifndef FUNCTION_SEAL_H
#define FUNCTION_SEAL_H

#include <vector>
#include "seal/seal.h"

using namespace seal;
using namespace std;

Ciphertext sgn(vector<double> input, vector<double> poly, int scale, int d, string mode, CKKSEncoder& encoder, Encryptor& enc, Evaluator& eva, SEALContext& context, RelinKeys& rlk, Decryptor& dec);
Ciphertext abs_seal(vector<double> input, vector<double> poly, int scale, int d, CKKSEncoder& encoder, Encryptor& enc, Evaluator& eva, SEALContext& context, RelinKeys& rlk, Decryptor& dec);
Ciphertext minMax_seal(vector<double> input, vector<double> poly, int scale, int d, string mode, CKKSEncoder& encoder, Encryptor& enc, Evaluator& eva, SEALContext& context, RelinKeys& rlk, Decryptor& dec);

#endif // FUNCTION_SEAL_H