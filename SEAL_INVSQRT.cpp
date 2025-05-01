#include "header/SEAL_VS.h"

/*
	Todo
	1. newton method(plain, ctxt)
	2. goldschmidt algorithm(plain, ctxt)
	3. good initial guess algorithm(plain)
	4. good initial guess algorithm(ctxt)
	5. depth, time consumption comparison
*/
int main() {
	double x = 2;
	cout << "Initial X : " << x << endl;

	string printmode = "debug";
	double delta = 1e-5;
	double err = 1e-4;
	int iter = 9;
	
	// k1, k2 계산
	auto [k1, k2] = find_bounds(iter, delta, err, printmode);
	debug_print(format("k1: {}", k1), printmode);
	debug_print(format("k2: {}", k2), printmode);

	debug_print("---------------------------------", printmode);

	// L2(x2) -> P -> L1(x1) 계산
	double a = 1e-4;
	double b = 1e3; 

	double x2 = solve_eq5(k1, k2, b, 400);
	vector<double> L2 = tangent_coeff(k1, k2, x2);
	debug_print(format("Tangent point x2:\t{}", x2), printmode);

	double pivot = solve_eq5_for_x(k1, k2, x2, 1.0);
	debug_print(format("Pivot point P:\t\t{}", pivot), printmode);

	double x1 = solve_eq5(k1, k2, pivot, 1.0);
	vector<double> L1 = tangent_coeff(k1, k2, x1);
	debug_print(format("Tangent point x1:\t{}", x1), printmode);

	debug_print("---------------------------------", printmode);

	debug_print("Tangent function L1:\t", printmode);
	printVector(L1, true);
	debug_print("Tangent function L2:\t", printmode);
	printVector(L2, true);
	
	debug_print("---------------------------------", printmode);

	// approximate invsqrt(x)	
	double real_value = 1 / sqrt(x);
	for (int i = 0; i < 100; i++) {
		//double x0 = sample_data(x1, 3)[0]; // sample x0 in [x1, x2]
		//double y0 = compute_h(pivot, x0, dg, df, L1, L2);
		//debug_print(format("x0:\t{}", x0), printmode);
		double y0 = sample_data(k1 / sqrt(x), k2 / sqrt(x))[0];
		debug_print(format("y0:\t{}", y0), printmode);

		double yi = y0;
		for (int i = 0; i < iter; i++) {
			yi = newton(x, yi);
			debug_print(format("y{}:\t{}", (i + 1), yi), printmode);
			if (abs(yi - real_value) <= err)
				break;
		}
		debug_print("---------------------------------", printmode);
	}

	return 0;
}
//int main()
//{
//	// param settings
//	int num = 5; // samples
//	int n = 1;
//	int nf = 5; // f_n
//	int ng = 5; // g_n
//	int alpha = 7;
//	double epsilon = pow(2.0, -alpha);
//	int pmd = 32768;
//	int big_moduli = 60;
//	int small_moduli = 60;
//	double scale = pow(2.0, 25);
//	int depth = 11;
//	ckks_build ckks = ckks_build(n, depth, big_moduli, small_moduli, scale, pmd);
//
//	// Sample data in [-1, -e] U [e, 1].
//	vector<double> raw_inputs = sample_data(-1.0, 1.0, epsilon, num);
//	cout << "sampled x\n\t";
//	printVector(raw_inputs, false, 6);
//	cout << "---" << endl;
//
//	Ciphertext x = ckks.encrypt(raw_inputs);
//
//	// compute f_n coeffs.
//	vector<double> polyF = computeF(n);
//	// use g_n based on: Cheon, Jung Hee, Dongwoo Kim, and Duhyeong Kim. "Efficient homomorphic comparison methods with optimal complexity."
//	vector<double> polyG = { 0, 2126 / pow(2.0, 10), 0,  -1359 / pow(2.0, 10) }; //g_1
//	cout << "f_" << n << ":\t";
//	printVector(polyF, true);
//	cout << "g_1" << ":\t";
//	printVector(polyG, true);
//	cout << "---" << endl;
//
//	auto time1 = cur_time(); //start time
//
//	vector<double> realValue = raw_inputs;
//	vector<double> fnDec;
//	// evaluate
//	Ciphertext y = x;
//	cout << "Initial Level: " << y.coeff_modulus_size() << endl;
//	cout << "Initial Scale: " << y.scale() << endl;
//	cout << "------------------------------------" << endl;
//
//	// g_1^ng
//	for (int i = 1; i <= ng; i++) {
//		ckks.evaluate_function_tripleScale_v2(polyG, y, y);
//		printStep(realValue, polyG, fnDec, raw_inputs, y, ckks, "g", i);
//	}
//
//	// f_1^nf
//	for (int i = 1; i <= nf; i++) {
//		ckks.evaluate_function_tripleScale_v2(polyF, y, y);
//		printStep(realValue, polyF, fnDec, raw_inputs, y, ckks, "f", i);
//	}
//
//	auto time2 = cur_time(); // end time
//	calculate_time(time1, time2);
//
//	return 0;
//}
