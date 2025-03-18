#include "header/SEAL_VS.h"

int main()
{
	// param settings
	string mode = "debug";
	int num = 5; // samples
	int n = 1; // f_n
	int alpha = 8;
	double epsilon = pow(2.0, -alpha);
	int big_moduli = 30;
	int small_moduli = 23;
	double scale = pow(2.0, small_moduli);
	//size_t pmd = 32768;
	size_t pmd = 16384;
	double c_n;
	int d, depth;
	int pre = 10;

	//choose function
	int function_mode;
	cout << "1.sgn(x), 2.abs(x), 3.min/max(x)" << endl;
	cin >> function_mode;
	cout << "---" << endl;

	// calculate d, total depth
	c_n = ((2 * n + 1) / pow(4.0, n)) * (factorial(2 * n, n) / factorial(n));
	double a = alpha / log2(c_n);
	double b = log2(alpha - 1) / log2(n + 1);
	d = a+b; // O(1)=1
	switch (function_mode)
	{
	case 1: // sgn
		depth = d * (2 * n);
		break;
	case 2: //abs
		depth = d * (2 * n) + 1;
		break;
	case 3: //max
		depth = d * (2 * n) + 2;
		break;
	default:
		return 0;
	}
	if (mode == "debug") {
		cout << "c_n = " << c_n << endl;
		cout << "d = " << a << "+" << b << "= " << d << endl;
		cout << "depth = " << depth << endl;
	}
	d = 4;
	depth = 11;

	ckks_build ckks = ckks_build(mode, n, depth, big_moduli, small_moduli, scale, pmd);

	//Sample data in [-1, -e] U [e, 1] and encrypt.
	vector<double> raw_inputs = sample_data(-1.0, 1.0, epsilon, num);
	cout << "sampled x\n\t";
	printVector(raw_inputs, false, 6);
	vector<double> realValue;
	realValue.resize(raw_inputs.size());
	Ciphertext result;
	Ciphertext x = ckks.encrypt(raw_inputs);
	cout << "---" << endl;

	//compute f_n coeffs.
	vector<double> poly = computeF(n);
	cout << "f_" << n << ":\t";
	printVector(poly, true);
	cout << "---" << endl;

	time_point<high_resolution_clock> start = cur_time();
	time_point<high_resolution_clock> end;
	switch (function_mode) 
	{
	case 1: // sgn(x)
		for (int i = 0; i < raw_inputs.size(); i++) {
			realValue[i] = polypolyEvaluate(poly, raw_inputs[i], d);
		}
		result = sgn_seal(mode, x, poly, raw_inputs, d, pre, ckks);
		calculate_time(start, cur_time());
		/*printResult(result, realValue, raw_inputs.size(), ckks, pre);*/
		cout << "---" << endl;
		break;

	case 2: // abs(x)
		for (int i = 0; i < raw_inputs.size(); i++) {
			realValue[i] = calAbs(raw_inputs[i], n, d);
		}
		result = abs_seal(mode, x, poly, raw_inputs, d, pre, ckks);
		calculate_time(start, cur_time());
		printResult(result, realValue, raw_inputs.size(), ckks, pre);
		cout << "---" << endl;
		break;

	case 3: // max(x)
		double maxValue;
		Ciphertext result;
		if (raw_inputs.size() == 2) {
			maxValue = raw_inputs[0] > raw_inputs[1] ? raw_inputs[0] : raw_inputs[1];
			result = max_seal(mode, x, poly, raw_inputs, d, pre, ckks);
			calculate_time(start, cur_time());
		}
		break;
	}

	//복호화
	vector<double> plain_outputs = ckks.decode_ctxt(result);
	plain_outputs.resize(realValue.size());
	cout << "Decryption Completed." << endl << "---" << endl;

	//결과 출력
	cout << "암호화 근사 결과" << endl;
	printVector(plain_outputs, false, 10);
	cout << "비암호화 근사 결과" << endl;
	printVector(realValue, false, 10);
	cout << "오차" << endl;
	vector<double> error;
	for (int i = 0; i < plain_outputs.size(); i++) {
		error.push_back(abs(plain_outputs[i] - realValue[i]));
	}
	printVector(error, false, 10);
	cout << endl;

	

	return 0;
}

//test double Scale
//int main()
//{
//	// param settings
//	string mode = "d";
//	int num = 5; // samples
//	int n = 1; // f_n
//	int alpha = 8;
//	double epsilon = pow(2.0, -alpha);
//	int big_moduli = 60;
//	int small_moduli = 40;
//	double scale = pow(2.0, small_moduli / 2);
//	size_t pmd = 32768;
//	int depth = 3;
//
//	ckks_build ckks = ckks_build(mode, n, depth, big_moduli, small_moduli, scale, pmd);
//
//	//Sample data in [-1, -e] U [e, 1].
//	vector<double> raw_inputs = sample_data(-1.0, 1.0, epsilon, num);
//	cout << "sampled x\n\t";
//	printVector(raw_inputs, false, 6);
//	cout << "---" << endl;
//
//	vector<double> realValue;
//	realValue.resize(raw_inputs.size());
//	Ciphertext result;
//	Ciphertext x = ckks.encrypt(raw_inputs);
//
//	//x = input samples, c = 1.0
//	Plaintext c = ckks.encode(1.0);
//	Ciphertext term0, term1, term2;
//
//	//evaluate (c^2)x
//	ckks.mul_plain(c, x, term0);  // s , L-1
//	c = ckks.encode(1.0, term0);
//	ckks.mul_plain(c, term0, term0);  // s, L-2
//
//	auto result_vec1 = ckks.decode_ctxt(term0);
//	result_vec1.resize(raw_inputs.size());
//	cout << "1.0x^2\ndecryption result:\t";
//	printVector(result_vec1, false, 6);
//
//	//evaluate x^3
//	ckks.mul_cipher(x, x, x, term1); // s, L-1
//	c = ckks.encode(1.0, term1);
//	ckks.mul_plain(c, term1, term1); // s, L-2
//
//	auto result_vec2 = ckks.decode_ctxt(term1);
//	result_vec2.resize(raw_inputs.size());
//	cout << "x^3\ndecryption result:\t";
//	printVector(result_vec2, false, 6);
//
//	//evaluate 1.0x + x^3
//	ckks.add_cipher(term0, term1, term2);
//
//	auto result_vec3 = ckks.decode_ctxt(term1);
//	result_vec3.resize(raw_inputs.size());
//	cout << "1.0x+x^3\ndecryption result:\t";
//	printVector(result_vec3, false, 6);
//
//	return 0;
//}
