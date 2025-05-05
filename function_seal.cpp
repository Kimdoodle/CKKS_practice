#include "header/SEAL_VS.h"
using namespace std;

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

// newton method
Ciphertext newton_seal(double x, vector<double> x0, int iter, string printmode, ckks_build& ckks)
{
	vector<double> y0 = x0;
	Ciphertext y = ckks.encrypt(x0);

	//iter newton method.
	Ciphertext term;
	Plaintext cipher_x, plain_3, plain_05, dummy;

	for (int i = 0; i < iter; i++) {
		debug_print(format("Iteration {}", i + 1), printmode);

		ckks.mult(y, y, term, false); // y^2

		cipher_x = ckks.encode((-1) * x, y);
		ckks.mult(cipher_x, term, false); //-xy^2

		plain_3 = ckks.encode(3.0, term);
		ckks.add(plain_3, term);// 3+(-xy^2)

		plain_05 = ckks.encode(0.5, term, pow(2.0, 10));
		ckks.mult(plain_05, term, true); // 0.5 * (3-xy^2)

		ckks.mult(y, term, term, false);// 0.5 * y * (3-xy^2)
		dummy = ckks.encode(1.0, term, pow(2.0, 35));
		ckks.mult(dummy, term, true);

		y = term;

		int iteration_scale = int(log2(y.scale()));
		if (iteration_scale != ckks.get_scale())
		{
			cout << "SCALE ERROR!!!" << endl;
			cout << format("iteration scale: {}", iteration_scale) << endl;
			cout << format("original scale: {}", ckks.get_scale()) << endl;
			return term;
		}

		//error
		if (printmode == "debug")
		{
			int size = x0.size();
			vector<double> res_ctxt = ckks.decode_ctxt(term);
			vector<double> res_error;
			res_ctxt.resize(x0.size());
			for (int i = 0; i < size; i++)
			{
				res_error.push_back(abs((1/sqrt(x)) - res_ctxt[i]));
			}
			cout << "Decryption Result:" << endl;
			printVector(res_ctxt, false);
			cout << "Error:" << endl;
			printVector_10eform(res_error, false);
			
		}

		cout << format("Remain Level: {}", term.coeff_modulus_size()) << endl;
		cout << "------------------------------" << endl;
	}
	
	return term;
}

Ciphertext goldschmidt_seal(double x, vector<double> x0, int iter, string printmode, ckks_build& ckks)
{
	vector<double> y0 = x0;
	Ciphertext y = ckks.encrypt(x0);

	//iter goldschmidt method.
	Plaintext cipher_x = ckks.encode((-1) * x, y);
	Ciphertext h;
	Plaintext plain_3, plain_05, dummy;
	
	//set initial y, g
	ckks.mult(y, y, h, false);
	ckks.mult(cipher_x, h, false);
	plain_3 = ckks.encode(3.0, h);
	ckks.add(plain_3, h);
	plain_05 = ckks.encode(0.5, h, pow(2.0, 10));
	ckks.mult(plain_05, h, true); //h0

	for (int i = 0; i < iter; i++) {
		debug_print(format("Iteration {}", i + 1), printmode);
		
		ckks.mult(y, h, y, false); // y_(n+1) = y_n * h_n
		dummy = ckks.encode(1.0, y, pow(2.0, 35));
		ckks.mult(dummy, y, true);

		ckks.mult(h, h, h, false); // h_(n+1) = h_n * h_n
		dummy = ckks.encode(1.0, h, pow(2.0, 35));
		ckks.mult(dummy, h, true);

		int iteration_scale = int(log2(y.scale()));
		if (iteration_scale != ckks.get_scale())
		{
			cout << "SCALE ERROR!!!" << endl;
			return y;
		}

		//error
		if (printmode == "debug")
		{
			int size = x0.size();
			vector<double> res_ctxt = ckks.decode_ctxt(y);
			vector<double> res_error;
			res_ctxt.resize(x0.size());
			for (int i = 0; i < size; i++)
			{
				res_error.push_back(abs((1 / sqrt(x)) - res_ctxt[i]));
			}
			cout << "Decryption Result:" << endl;
			printVector(res_ctxt, false);
			cout << "Error:" << endl;
			printVector_10eform(res_error, false);

		}

		cout << format("Remain Level: {}", y.coeff_modulus_size()) << endl;
		cout << "------------------------------" << endl;
	}

	return y;
}