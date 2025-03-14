#ifndef COMPARE_H
#define COMPARE_H

#include "SEAL_VS.h"

double calC(int n);
std::vector<double> computeF(int n);
std::vector<double> computeH(int n);
double signFunction(double a, int d);
double newComp(double a, double b, int n, int d);
double calAbs(double a, int n, int d);
double calMin(double a, double b, int n, int d);
double calMax(double a, double b, int n, int d);

#endif // COMPARE_H
