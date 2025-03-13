#ifndef POLY_EVAL_H
#define POLY_EVAL_H

#include <vector>
#include "seal/seal.h"

using namespace seal;
using namespace std;

Ciphertext exp_x(Ciphertext x, int d, Evaluator& eva, SEALContext& context, RelinKeys& rlk, double scale);
Ciphertext mult_x_plain(Plaintext coeff, Ciphertext x, double scale, Evaluator& eva);
Ciphertext cal_x1_x2(char mode, Ciphertext x1, Ciphertext x2, const SEALContext& context, Evaluator& eva, RelinKeys& rlk, double scale);

#endif // POLY_EVAL_H
