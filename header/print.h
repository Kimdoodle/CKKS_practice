#ifndef PRINT_H
#define PRINT_H

#include "SEAL_VS.h"

void printVector(std::vector<double> coeffs, bool asFunction);
void print_parameters(const seal::SEALContext& context);

#endif // PRINT_H
