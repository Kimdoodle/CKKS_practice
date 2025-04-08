#include "header/SEAL_VS.h"

int main()
{
	// param settings
	int num = 5; // samples
	int n = 1;
	int nf = 5; // f_n
	int ng = 5; // g_n
	int alpha = 7;
	double epsilon = pow(2.0, -alpha);
	int pmd = 32768;
	int big_moduli = 60;
	int small_moduli = 60;
	double scale = pow(2.0, 25);
	int depth = 11;
	ckks_build ckks = ckks_build(n, depth, big_moduli, small_moduli, scale, pmd);

	// Sample data in [-1, -e] U [e, 1].
	vector<double> raw_inputs = sample_data(-1.0, 1.0, epsilon, num);
	cout << "sampled x\n\t";
	printVector(raw_inputs, false, 6);
	cout << "---" << endl;

	Ciphertext x = ckks.encrypt(raw_inputs);

	// compute f_n coeffs.
	vector<double> polyF = computeF(n);
	// use g_n based on: Cheon, Jung Hee, Dongwoo Kim, and Duhyeong Kim. "Efficient homomorphic comparison methods with optimal complexity."
	vector<double> polyG = { 0, 2126 / pow(2.0, 10), 0,  -1359 / pow(2.0, 10) }; //g_1
	cout << "f_" << n << ":\t";
	printVector(polyF, true);
	cout << "g_1" << ":\t";
	printVector(polyG, true);
	cout << "---" << endl;

	auto time1 = cur_time(); //start time

	vector<double> realValue = raw_inputs;
	vector<double> fnDec;
	// evaluate
	Ciphertext y = x;
	cout << "Initial Level: " << y.coeff_modulus_size() << endl;
	cout << "Initial Scale: " << y.scale() << endl;
	cout << "------------------------------------" << endl;

	// g_1^ng
	for (int i = 1; i <= ng; i++) {
		ckks.evaluate_function_tripleScale_v2(polyG, y, y);
		printStep(realValue, polyG, fnDec, raw_inputs, y, ckks, "g", i);
	}

	// f_1^nf
	for (int i = 1; i <= nf; i++) {
		ckks.evaluate_function_tripleScale_v2(polyF, y, y);
		printStep(realValue, polyF, fnDec, raw_inputs, y, ckks, "f", i);
	}

	auto time2 = cur_time(); // end time
	calculate_time(time1, time2);

	return 0;
}
