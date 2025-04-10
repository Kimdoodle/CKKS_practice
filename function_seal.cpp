#include "header/SEAL_MM.h"

//evaluate function
Ciphertext evaluate_function(vector<double>& poly, Ciphertext& x, int size, string scaleMode, ckks_build& ckks)
{
	Ciphertext result = ckks.encrypt(0.0);
	for (int i = 0; i < poly.size(); i++) {
		if (poly[i] == 0.0) continue;
		Ciphertext term1 = x;
		Plaintext coeff = ckks.encode(poly[i], x); // a
		ckks.mult(coeff, term1); //ax	
		Ciphertext squareX = ckks.exp(x, i); // x^(i-1)
		ckks.mult(term1, squareX); // ax^i
		ckks.add(result, squareX); // y += ax^i
	}
	return result;
}

// sign function
Ciphertext sgn_seal(string mode, string scaleMode, Ciphertext& x, vector<double> poly, vector<double>& input, int d, int pre, ckks_build& ckks)
{
	Ciphertext y = x;
	vector<double> real_result = input;
	if (mode == "debug") {
		cout << "\tRemaining Level: " << y.coeff_modulus_size() << endl;
	}
	for (int iter = 0; iter < d; iter++) {
		y = evaluate_function(poly, y, input.size(), scaleMode, ckks); //y=f(x)
	}
	return y;
}

// abs function
Ciphertext abs_seal(string mode, string scaleMode, Ciphertext& x, vector<double>& poly, vector<double>& input, int d, int pre, ckks_build& ckks)
{
	Ciphertext result;
	Ciphertext sgnx = sgn_seal(mode, scaleMode, x, poly, input, d, pre, ckks);
	ckks.mult(x, sgnx, result);
	return result;
}

// max function
Ciphertext max_seal(string mode, string scaleMode, Ciphertext& x, vector<double>& poly, vector<double>& input, int d, int pre, ckks_build& ckks)
{
	Ciphertext ssum_abs_half, mmin_abs_half;
	Ciphertext ctxt, abs_result, half;
	Ciphertext result;
	
	// |a+b|/2
	double ssum = input[0] + input[1];
	ctxt = ckks.encrypt(ssum);
	abs_result = abs_seal(mode, scaleMode, x, poly, input, d, pre, ckks);
	half = ckks.encrypt(0.5);
	ckks.mult(half, abs_result, ssum_abs_half);
	
	// |a-b|/2
	double mmin = input[0] - input[1];
	ctxt = ckks.encrypt(mmin);
	abs_result = abs_seal(mode, scaleMode, x, poly, input, d, pre, ckks);
	ckks.mult(half, abs_result, mmin_abs_half);

	// calculate
	ckks.add(ssum_abs_half, mmin_abs_half, result);
	return result;
}