#include "header/SEAL_VS.h"

//evaluate function
Ciphertext evaluate_function(vector<double>& poly, Ciphertext& x, int size, ckks_build& ckks)
{
	Ciphertext result = ckks.encrypt(0.0);

	for (int i = 0; i < poly.size(); i++) {
		if (poly[i] == 0.0) continue;
		Ciphertext squareX = ckks.exp(x, i); // x^i
		Plaintext coeff = ckks.encode(poly[i], squareX); // a
		ckks.mult(coeff, squareX); // ax^i
		ckks.add(result, squareX); // y += ax^i
	}
	return result;
}

// sign function
Ciphertext sgn(string mode, Ciphertext& x, vector<double>& poly, vector<double>& input, int d, ckks_build& ckks)
{
	Ciphertext y = x;
	vector<double> real_result = input;
	if (mode == "debug") {
		cout << "\tRemaining Level: " << y.coeff_modulus_size() << endl;
	}
	for (int iter = 0; iter < d; iter++) {

		y = evaluate_function(poly, y, input.size(), ckks); //y=f(x)

		if (mode == "debug") {
			cout << "d = " << iter+1 << endl;
			//Remaining Levels
			cout << "\tRemaining Level: " << y.coeff_modulus_size() << endl;
			//복호화
			cout << "\tdecrypt value\n\t";
			vector<double> result = ckks.decode_ctxt(y);
			result.resize(input.size());
			printVector(result, false);
			//실제값
			real_result = polypolyEvaluate(poly, real_result);
			cout << "\treal value\n\t";
			printVector(real_result, false);
		}
	}
	return y;
}

// abs function
// function calculates all abs value in input.
//Ciphertext abs_seal(Ciphertext& x, vector<double>& poly, vector<double>& input, int d, ckks_build& ckks)
//{
//	Ciphertext result;
//	Ciphertext sgnx = sgn("debug", x, poly, input, d, ckks);
//	ckks.mult(x, sgnx, result);
//
//	return result;
//}

// max function
//Ciphertext max_seal(vector<double>& input, vector<double>& poly, int d, ckks_build& ckks)
//{
//	Ciphertext ssum_abs_half, mmin_abs_half;
//	Ciphertext ctxt, abs, half;
//	Ciphertext result;
//	
//	// |a+b|/2
//	double ssum = input[0] + input[1];
//	ctxt = ckks.encrypt(ssum);
//	abs = abs_seal(ctxt, poly, d, ckks);
//	half = ckks.encrypt(0.5);
//	ckks.mult(half, abs, ssum_abs_half);
//	
//	// |a-b|/2
//	double mmin = input[0] - input[1];
//	ctxt = ckks.encrypt(mmin);
//	abs = abs_seal(ctxt, poly, d, ckks);
//	ckks.mult(half, abs, mmin_abs_half);
//
//	// calculate
//	ckks.add(ssum_abs_half, mmin_abs_half, result);
//	return result;
//}