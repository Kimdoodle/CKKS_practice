#include "header/SEAL_VS.h"

//evaluate function
Ciphertext evaluate_function(vector<double>& poly, Ciphertext& x, ckks_build& ckks) 
{
	Ciphertext result = ckks.encrypt(0.0);

	for (int i = 0; i < poly.size(); i++) {
		if (poly[i] == 0.0) continue;
		Ciphertext squareX = ckks.exp(x, i);
		Plaintext coeff = ckks.encode(poly[i], squareX);
		ckks.mult(coeff, squareX);
		ckks.add(result, squareX);
	}
	return result;
}
// sign function
Ciphertext sgn(string mode, Ciphertext& x, vector<double>& poly, int d, ckks_build& ckks)
{
	Ciphertext y = x;
	for (int iter = 0; iter < d; iter++) {
		y = evaluate_function(poly, y, ckks);
	}
	return y;
}

// abs function
// function calculates all abs value in input.
Ciphertext abs_seal(Ciphertext& input, vector<double>& poly, int d, ckks_build& ckks)
{
	Ciphertext result;
	Ciphertext x = input;
	Ciphertext sgnx = sgn("debug", input, poly, d, ckks);
	ckks.mult(x, sgnx, result);

	return result;
}

// max function
Ciphertext max_seal(vector<double>& input, vector<double>& poly, int d, ckks_build& ckks)
{
	Ciphertext ssum_abs_half, mmin_abs_half;
	Ciphertext ctxt, abs, half;
	Ciphertext result;
	
	// |a+b|/2
	double ssum = input[0] + input[1];
	ctxt = ckks.encrypt(ssum);
	abs = abs_seal(ctxt, poly, d, ckks);
	half = ckks.encrypt(0.5);
	ckks.mult(half, abs, ssum_abs_half);
	
	// |a-b|/2
	double mmin = input[0] - input[1];
	ctxt = ckks.encrypt(mmin);
	abs = abs_seal(ctxt, poly, d, ckks);
	ckks.mult(half, abs, mmin_abs_half);

	// calculate
	ckks.add(ssum_abs_half, mmin_abs_half, result);
	return result;
}