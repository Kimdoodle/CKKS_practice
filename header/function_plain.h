#pragma once
#include "SEAL_VS.h"

//double calC(int n);
//vector<double> computeF(int n);
//vector<double> computeH(int n);
//vector<double> computeG(int n, double tau, double pre, double a, double b);
//double signFunction(double a, int d);
//double newComp(double a, double b, int n, int d);
//double calAbs(double a, int n, int d);
//double calMin(double a, double b, int n, int d);
//double calMax(double a, double b, int n, int d);

//invsqrt functions
double newton(double x, double y);
double iter_newton(double x, double y, int iter);
double calculate_k1(double low, double high, int iter, double delta, double err, string printmode);
double calculate_k2(double low, double high, int iter, double delta, double err, string printmode);
pair<double, double> find_bounds(int iter, double delta, double err, const string& mode);

double solve_eq5(double k1, double k2, double x, double x0, int max_iter = 100, double tol = 1e-10);
double solve_eq5_for_x(double k1, double k2, double x2_fixed, double x0, int max_iter = 100, double tol = 1e-10);
vector<double> tangent_coeff(double k1, double k2, double x);

double approx_sign(double x, int dg, int df);
double approx_comp(double a, double b, int dg, int df);
double compute_h(double p, double x, int dg, int df, vector<double> L1, vector<double> L2);

int newton_algorithm(double x, double x0, double err, int iter, string printmode);
int goldschmidt_algorithm(double x, double x0, double err, int iter, string printmode);
void newton_goodGuess(double x, double a, double b, double delta, double err, int iter, string printmode);