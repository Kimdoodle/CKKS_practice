#include "header/SEAL_VS.h"

//int main()
//{
//	// param settings
//	string mode = "debug";
//	int num = 5; // samples
//	int n = 1; // f_n
//	int alpha = 8;
//	double epsilon = pow(2.0, -alpha);
//	int big_moduli = 30;
//	int small_moduli = 23;
//	double scale = pow(2.0, small_moduli);
//	//size_t pmd = 32768;
//	size_t pmd = 16384;
//	double c_n;
//	int d, depth;
//	int pre = 10;
//
//	//choose function
//	int function_mode;
//	cout << "1.sgn(x), 2.abs(x), 3.min/max(x)" << endl;
//	cin >> function_mode;
//	cout << "---" << endl;
//
//	// calculate d, total depth
//	c_n = ((2 * n + 1) / pow(4.0, n)) * (factorial(2 * n, n) / factorial(n));
//	double a = alpha / log2(c_n);
//	double b = log2(alpha - 1) / log2(n + 1);
//	d = a+b; // O(1)=1
//	switch (function_mode)
//	{
//	case 1: // sgn
//		depth = d * (2 * n);
//		break;
//	case 2: //abs
//		depth = d * (2 * n) + 1;
//		break;
//	case 3: //max
//		depth = d * (2 * n) + 2;
//		break;
//	default:
//		return 0;
//	}
//	if (mode == "debug") {
//		cout << "c_n = " << c_n << endl;
//		cout << "d = " << a << "+" << b << "= " << d << endl;
//		cout << "depth = " << depth << endl;
//	}
//	d = 4;
//	depth = 11;
//
//	ckks_build ckks = ckks_build(mode, n, depth, big_moduli, small_moduli, scale, pmd);
//
//	//Sample data in [-1, -e] U [e, 1] and encrypt.
//	vector<double> raw_inputs = sample_data(-1.0, 1.0, epsilon, num);
//	cout << "sampled x\n\t";
//	printVector(raw_inputs, false, 6);
//	vector<double> realValue;
//	realValue.resize(raw_inputs.size());
//	Ciphertext result;
//	Ciphertext x = ckks.encrypt(raw_inputs);
//	cout << "---" << endl;
//
//	//compute f_n coeffs.
//	vector<double> poly = computeF(n);
//	cout << "f_" << n << ":\t";
//	printVector(poly, true);
//	cout << "---" << endl;
//
//	time_point<high_resolution_clock> start = cur_time();
//	time_point<high_resolution_clock> end;
//	switch (function_mode) 
//	{
//	case 1: // sgn(x)
//		for (int i = 0; i < raw_inputs.size(); i++) {
//			realValue[i] = polypolyEvaluate(poly, raw_inputs[i], d);
//		}
//		result = sgn_seal(mode, x, poly, raw_inputs, d, pre, ckks);
//		calculate_time(start, cur_time());
//		/*printResult(result, realValue, raw_inputs.size(), ckks, pre);*/
//		cout << "---" << endl;
//		break;
//
//	case 2: // abs(x)
//		for (int i = 0; i < raw_inputs.size(); i++) {
//			realValue[i] = calAbs(raw_inputs[i], n, d);
//		}
//		result = abs_seal(mode, x, poly, raw_inputs, d, pre, ckks);
//		calculate_time(start, cur_time());
//		printResult(result, realValue, raw_inputs.size(), ckks, pre);
//		cout << "---" << endl;
//		break;
//
//	case 3: // max(x)
//		double maxValue;
//		Ciphertext result;
//		if (raw_inputs.size() == 2) {
//			maxValue = raw_inputs[0] > raw_inputs[1] ? raw_inputs[0] : raw_inputs[1];
//			result = max_seal(mode, x, poly, raw_inputs, d, pre, ckks);
//			calculate_time(start, cur_time());
//		}
//		break;
//	}
//
//	//복호화
//	vector<double> plain_outputs = ckks.decode_ctxt(result);
//	plain_outputs.resize(realValue.size());
//	cout << "Decryption Completed." << endl << "---" << endl;
//
//	//결과 출력
//	cout << "암호화 근사 결과" << endl;
//	printVector(plain_outputs, false, 10);
//	cout << "비암호화 근사 결과" << endl;
//	printVector(realValue, false, 10);
//	cout << "오차" << endl;
//	vector<double> error;
//	for (int i = 0; i < plain_outputs.size(); i++) {
//		error.push_back(abs(plain_outputs[i] - realValue[i]));
//	}
//	printVector(error, false, 10);
//	cout << endl;
//
//	
//
//	return 0;
//}

// Double Scale
int main()
{
	// param settings
	string mode = "d";
	int num = 5; // samples
	int n = 1;
	int nf = 3; // f_n
	int ng = 5; // g_n
	int alpha = 5;
	double epsilon = pow(2.0, -alpha);
	int big_moduli = 60;
	int small_moduli = 40;
	double scale = pow(2.0, small_moduli / 2);
	size_t pmd = 32768;
	int depth =17;
	// 120+21=141/3=47     59 41 41
	ckks_build ckks = ckks_build(mode, nf, depth, big_moduli, small_moduli, scale, pmd);

	// Sample data in [-1, -e] U [e, 1].
	vector<double> raw_inputs = sample_data(-1.0, 1.0, epsilon, num);
	cout << "sampled x\n\t";
	printVector(raw_inputs, false, 6);
	cout << "---" << endl;

	Ciphertext x = ckks.encrypt(raw_inputs);

	// compute f_n coeffs.
	vector<double> polyF = computeF(n);
	vector<double> polyG = { 0, 2126 / pow(2.0, 10), 0,  -1359 / pow(2.0, 10) }; //g_1
	cout << "f_" << n << ":\t";
	printVector(polyF, true);
	cout << "g_1" << ":\t";
	printVector(polyG, true);
	cout << "---" << endl;
	

	vector<double> realValue = raw_inputs;
	vector<double> fnDec;
	// evaluate
	Ciphertext y = x;
	cout << "Remaining Levels: " << y.coeff_modulus_size() << endl;
	// g_1^ng
	for (int i = 1; i <= ng; i++) {
		ckks.temp_d3_doubleScale(polyG, y, y);
		realValue = polypolyEvaluate(polyG, realValue);

		fnDec = ckks.decode_ctxt(y);
		fnDec.resize(raw_inputs.size());
		cout << "g^" << i << " 실제 계산값:\n\t";
		printVector(realValue, false);
		cout << "g^" << i << " 복호화값:\n\t";
		printVector(fnDec, false);
		cout << "Remaining Levels: " << y.coeff_modulus_size() << endl;
		cout << "---" << endl;
	}

	// f_1^nf
	for (int i = 1; i <= nf; i++) {
		ckks.temp_d3_doubleScale(polyF, y, y);
		realValue = polypolyEvaluate(polyF, realValue);

		fnDec = ckks.decode_ctxt(y);
		fnDec.resize(raw_inputs.size());
		cout << "f^"<<i<<" 실제 계산값:\n\t";
		printVector(realValue, false);
		cout << "f^" << i << " 복호화값:\n\t";
		printVector(fnDec, false);
		cout << "Remaining Levels: " << y.coeff_modulus_size() << endl;
		cout << "---" << endl;
	}
	
	//// g_1 = ax + bx^3
	//coeff = ckks.encode(polyF[1], x);
	//Ciphertext term0, term1, result;
	//ckks.mul_plain_double(coeff, x, term0); // ax
	//coeff = ckks.encode(polyF[3], x);
	//ckks.mul_plain_double(coeff, x, term1);
	//ckks.mul_cipher_double(term1, x, x, term1); //bx^3
	//ckks.add_cipher(term0, term1, result); // ax + bx^3

	return 0;
}

//int main() {
//	int n = 1;
//	double tau = 1 / pow(2.0, 2);
//	double pre = pow(2.0, -10);
//	double a = pow(2.0, -8);
//	double b = 1;
//	
//	vector<double> result = computeG(n, tau, pre, a, b);
//	printVector(result, true, 6);
//}
