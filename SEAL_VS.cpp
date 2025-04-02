#include "header/SEAL_VS.h"

void gogo(vector<double>& realValue, vector<double>& poly, vector<double>& fnDec, vector<double> raw_inputs, Ciphertext& y, ckks_build& ckks, int i)
{
	realValue = polypolyEvaluate(poly, realValue);
	fnDec = ckks.decode_ctxt(y);
	fnDec.resize(raw_inputs.size());
	cout << "g^" << i << " 실제 계산값:\n\t";
	printVector(realValue, false, 6);
	cout << "g^" << i << " 복호화값:\n\t";
	printVector(fnDec, false, 6);
	cout << "Remaining Levels: " << y.coeff_modulus_size() << endl;
	cout << "---" << endl;
}


int main()
{
	// param settings
	string mode = "d";
	string scaleMode;
	int num = 5; // samples
	int n = 1;
	
	//printf("#######SCALEMODE: DOUBLE########\n");
	scaleMode = "double";
	int nf = 3; // f_n
	int ng = 5; // g_n
	int alpha = 5;
	double epsilon = pow(2.0, -alpha);
	size_t pmd = 32768;
	int big_moduli = 60;
	int small_moduli = 40;
	double scale = pow(2.0, small_moduli / 2);
	int depth = 17;

	//ckks_build ckks = ckks_build(mode, scaleMode, nf, depth, big_moduli, small_moduli, scale, pmd);

	// Sample data in [-1, -e] U [e, 1].
	//vector<double> raw_inputs = sample_data(-1.0, 1.0, epsilon, num);
	vector<double> raw_inputs;
	for (int i = 0; i < 5; i++)
		raw_inputs.push_back(0.0);
	//cout << "sampled x\n\t";
	//printVector(raw_inputs, false, 6);
	//cout << "---" << endl;

	Ciphertext x;
	//Ciphertext x = ckks.encrypt(raw_inputs);

	//// compute f_n coeffs.
	vector<double> polyF = computeF(n);
	vector<double> polyG = { 0, 2.0762, 0,  -1.3271 }; //g_1
	//vector<double> polyG = { 0, 2126 / pow(2.0, 10), 0,  -1359 / pow(2.0, 10) }; //g_1
	
	//cout << "f_" << n << ":\t";
	//printVector(polyF, true);
	//cout << "g_1" << ":\t";
	//printVector(polyG, true);
	//cout << "---" << endl;

	vector<double> realValue = raw_inputs;
	vector<double> fnDec;
	//// evaluate
	Ciphertext y = x;
	//cout << "Remaining Levels: " << y.coeff_modulus_size() << endl;
	//// g_1^ng
	//for (int i = 1; i <= ng; i++) {
	//	ckks.temp_d3_doubleScale(polyG, y, y);
	//	gogo(realValue, polyG, fnDec, raw_inputs, y, ckks, i);
	//}

	//// f_1^nf
	//for (int i = 1; i <= nf; i++) {
	//	ckks.temp_d3_doubleScale(polyF, y, y);
	//	gogo(realValue, polyF, fnDec, raw_inputs, y, ckks, i);
	//}


	printf("#######SCALEMODE: TRIPLE########\n");
	scaleMode = "triple";
	nf = 10; // f_n
	ng = 0; // g_n
	alpha = 7;
	epsilon = pow(2.0, -alpha);
	pmd = 32768;
	big_moduli = 60;
	small_moduli = 60;
	scale = pow(2.0, small_moduli / 3);
	depth = 11;
	ckks_build ckks = ckks_build(mode, scaleMode, n, depth, big_moduli, small_moduli, scale, pmd);

	// Sample data in [-1, -e] U [e, 1].
	//raw_inputs = sample_data(-1.0, 1.0, epsilon, num);
	cout << "sampled x\n\t";
	printVector(raw_inputs, false, 6);
	cout << "---" << endl;

	x = ckks.encrypt(raw_inputs);

	// compute f_n coeffs.
	//polyF = computeF(n);
	//polyG = { 0, 2126 / pow(2.0, 10), 0,  -1359 / pow(2.0, 10) }; //g_1
	cout << "f_" << n << ":\t";
	printVector(polyF, true);
	cout << "g_1" << ":\t";
	printVector(polyG, true);
	cout << "---" << endl;


	realValue = raw_inputs;
	fnDec;
	// evaluate
	y = x;
	cout << "Remaining Levels: " << y.coeff_modulus_size() << endl;
	// g_1^ng
	for (int i = 1; i <= ng; i++) {
		cout << "Scale: " << y.scale() << endl;
		ckks.temp_d3_tripleScale(polyG, y, y);
		//y.scale() = scale;
		gogo(realValue, polyG, fnDec, raw_inputs, y, ckks, i);
	}

	// f_1^nf
	for (int i = 1; i <= nf; i++) {
		cout << "Scale: " << y.scale() << endl;
		ckks.temp_d3_tripleScale(polyF, y, y);
		//y.scale() = scale;
		gogo(realValue, polyF, fnDec, raw_inputs, y, ckks, i);
	}

	
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
