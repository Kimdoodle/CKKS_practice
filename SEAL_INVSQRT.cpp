#include "header/SEAL_VS.h"
//void plain_compare();

/*
	Todo
	1. newton method(plain, ctxt)
	2. goldschmidt algorithm(plain, ctxt)
	3. good initial guess algorithm(plain) - clear
	4. good initial guess algorithm(ctxt)
	5. depth, time consumption comparison
*/
int main() 
{
	//params for newton method
	double err = 1e-4;
	double x = 2.0;
	double answer = 1 / sqrt(x);
	int x0_size = 10;
	vector<double> x0 = sample_data(0, sqrt(3) / sqrt(x), x0_size);	
	//vector<double> x0 = sample_data(0, 2, x0_size);	
	int iter = 5;
	string printmode = "debug";

	cout << "SAMPLED DATA" << endl;
	printVector(x0, false);
	cout << "#########################################################" << endl;

	//params for ckks
	int moduli = 60;
	double scale = pow(2.0, 25);
	size_t pmd = pow(2.0, 15);
	ckks_build ckks(moduli, scale, pmd);

	cout << "NEWTON METHOD" << endl;
	newton_seal(x, x0, iter, printmode, ckks);
	cout << "#########################################################" << endl;
	cout << "GOLDSCHMIDT METHOD(FHE)" << endl;
	goldschmidt_seal(x, x0, iter, printmode, ckks);
}


//void plain_compare() 
//{
//	double answer = 2.0;
//	double x0 = 1.0;
//	double err = 1e-4;
//	double iter = 9;
//	string printmode = "normal";
//
//	int suc_new;
//	int suc_gold;
//	int new_adv = 0;
//	int gold_adv = 0;
//	for (double a = 1e-5; a <= 1e4; a += 1e-5)
//	{
//		suc_new = newton_algorithm(answer, x0, err, iter, printmode);
//		suc_gold = goldschmidt_algorithm(answer, x0, err, iter, printmode);
//		if (suc_new == -1 && suc_gold != -1)
//		{
//			cout << format("X: {}", a) << endl;
//			cout << format("Newton failed.\t Goldschmidt success. Iter: {}", suc_gold) << endl;
//		}
//		else if (suc_new != -1 && suc_gold == -1)
//		{
//			cout << format("X: {}", a) << endl;
//			cout << format("Newton success. Iter: {}.\t Goldschmidt failed.", suc_new) << endl;
//		}
//		else if (suc_new > suc_gold)
//			new_adv++;
//		else if (suc_new < suc_gold)
//			gold_adv++;
//	}
//	cout << format("Newton advantage: {} times", new_adv) << endl;
//	cout << format("Goldschmidt advantage: {} times", gold_adv) << endl;
//}

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
