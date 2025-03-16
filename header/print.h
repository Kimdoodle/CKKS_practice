#ifndef PRINT_H
#define PRINT_H

#include "SEAL_VS.h"

void printVector(vector<double>& coeffs, bool asFunction, int pre = 4);
void print_parameters(const SEALContext& context);

#endif // PRINT_H
