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
void newton_seal(double x, vector<double> x0, int iter, string printmode, ckks_build& ckks)
{
	vector<double> res_ctxt;
	vector<double> plain_y = x0;
	Ciphertext y = ckks.encrypt(x0);

	//iter newton method.
	// y = 0.5y(3-xy^2)
	Ciphertext term;
	Plaintext encoded_x, plain_05, plain_10, plain_30;

	for (int i = 0; i < iter; i++) {
		debug_print(format("Iteration {}", i + 1), printmode);

		ckks.mult(y, y, term, false); // y^2

		encoded_x = ckks.encode((-1) * x, y);
		ckks.mult(encoded_x, term, term, false); //-xy^2

		plain_30 = ckks.encode(3.0, term);
		ckks.add(plain_30, term);// 3+(-xy^2)

		plain_10 = ckks.encode(1.0, term, pow(2, 10));
		ckks.mult(plain_10, term, true); // Rescale 3+(-xy^2)
		
		//if (printmode == "debug")
		//{
		//	vector<double> res_ctxt = ckks.decode_ctxt(term);
		//	vector<double> res_ptxt = multVectors(x0, x0);
		//	res_ptxt = multPlainPolynomial(res_ptxt, -x);
		//	res_ptxt = addScalar(res_ptxt, 3);
		//	res_ctxt.resize(res_ptxt.size());
		//	cout << "decyrption result of 3-xy^2:" << endl;
		//	printVector(res_ctxt, false);
		//	cout << "Actual 3-xy^2:" << endl;
		//	printVector(res_ptxt, false);
		//}

		int iteration_scale = int(log2(term.scale()));
		if (iteration_scale != ckks.get_scale())
		{
			cout << "SCALE ERROR!!!" << endl;
			cout << format("iteration scale: {}", iteration_scale) << endl;
			cout << format("original scale: {}", ckks.get_scale()) << endl;
		}
		//###################################################
		plain_05 = ckks.encode(0.5, term);
		ckks.mult(plain_05, term, term, false); // 0.5 * (3-xy^2)

		ckks.mult(y, term, term, false);// 0.5 * y * (3-xy^2)

		plain_10 = ckks.encode(1.0, term, pow(2, 10));
		ckks.mult(plain_10, term, true); // Rescale

		y = term;

		//error
		if (printmode == "debug")
		{
			int size = x0.size();
			res_ctxt = ckks.decode_ctxt(y);
			res_ctxt.resize(x0.size());
			plain_y = newton_algorithm(x, plain_y);
			vector<double> res_error = subVectors(res_ctxt, plain_y);
			vector<double> res_real_error = subScalar(res_ctxt, 1 / sqrt(x));
			res_real_error = absVectors(res_real_error);

			//cout << "Decryption Result:" << endl;
			//printVector(res_ctxt, false);
			//cout << "Actual calculation result:" << endl;
			//printVector(plain_y, false);
			//cout << "Error:" << endl;
			//printVector_10eform(res_error, false);
			//cout << "Error(with real invsqrt):" << endl;
			//printVector_10eform(res_real_error, false);
			cout << "Approximate success?(Plaintext).. \t";
			vector<double> res_err_plain = subScalar(plain_y, 1 / sqrt(x));
			res_err_plain = absVectors(res_err_plain);
			printVector_OX(res_err_plain, 1e-4);
			cout << "Approximate success?(Ciphertext)\t";
			printVector_OX(res_real_error, 1e-4);
		}

		//iteration_scale = int(log2(y.scale()));
		//if (iteration_scale != ckks.get_scale())
		//{
		//	cout << "SCALE ERROR!!!" << endl;
		//	cout << format("iteration scale: {}", iteration_scale) << endl;
		//	cout << format("original scale: {}", ckks.get_scale()) << endl;
		//}
		

		cout << format("Remain Level: {}", y.coeff_modulus_size()) << endl;
		cout << "------------------------------" << endl;

		//zzap-bootstrapping
		if (y.coeff_modulus_size() == 1 && i != iter-1)
		{
			cout << "Re-encrypting..." << endl;
			y = ckks.encrypt(res_ctxt);
			cout << "------------------------------" << endl;
		}
	}
}

//goldschmidt algorithm
void goldschmidt_seal(double x, vector<double> x0, int iter, string printmode, ckks_build& ckks)
{

	Ciphertext g, h, temp;
	Plaintext encoded_x, plain_30, plain_05, plain_10;
	vector<double> res_ctxt;
	vector<double> res_h(x0.size()), res_y, res_g(x0.size()), res_err, res_real_err;
	Ciphertext y = ckks.encrypt(x0);

	res_y = x0;
	int iteration_scale;

	//set initial g = xy^2
	ckks.mult(y, y, temp, false);
	encoded_x = ckks.encode(x, temp, pow(2, 35));
	ckks.mult(encoded_x, temp, g, true);
	//error
	if (printmode == "debug")
	{
		res_ctxt = ckks.decode_ctxt(g);
		res_ctxt.resize(x0.size());
		res_g = multVectors(res_y, res_y);
		res_g = multScalar(res_g, x);
		//cout << "Decryption result g0:" << endl;
		//printVector(res_ctxt, false);
		//cout << "Actual calculation of g0:" << endl;
		//printVector(res_g, false);
		//iteration_scale = int(log2(g.scale()));
		//if (iteration_scale != ckks.get_scale())
		//{
		//	cout << "G0 SCALE ERROR!!!" << endl;
		//	cout << format("iteration scale: {}", iteration_scale) << endl;
		//	cout << format("original scale: {}", ckks.get_scale()) << endl;
		//	return;
		//}
	}
	cout << "------------------------------" << endl;

	//iter goldschmidt method.
	for (int i = 0; i < iter; i++) {
		debug_print(format("Iteration {}", i + 1), printmode);

		//calculate with plain data
		if (printmode == "debug")
		{
			goldschmidt_algorithm(x, res_y, res_g, res_h);
		}

		//calculate h = (3-g)/2 = (g-3)/(-2)
		plain_30 = ckks.encode(-3.0, g);
		ckks.add(plain_30, g, temp);
		plain_05 = ckks.encode(-0.5, temp, pow(2, 60));
		ckks.mult(plain_05, temp, h, true);
		//error
		//if (printmode == "debug")
		//{
		//	res_ctxt = ckks.decode_ctxt(h);
		//	res_ctxt.resize(x0.size());
		//	cout << "Decryption Result of h:" << endl;
		//	printVector(res_ctxt, false);
		//	cout << "Actual calculation of h:" << endl;
		//	printVector(res_h, false);
		//	iteration_scale = int(log2(h.scale()));
		//	if (iteration_scale != ckks.get_scale())
		//	{
		//		cout << "H SCALE ERROR!!!" << endl;
		//		cout << format("iteration scale: {}", iteration_scale) << endl;
		//		cout << format("original scale: {}", ckks.get_scale()) << endl;
		//		return;
		//	}
		//	cout << "--------" << endl;
		//}

		//calculate g = gh^2
		ckks.mult(h, h, temp, false);
		ckks.mult(temp, g, temp, false);
		plain_10 = ckks.encode(1.0, temp, pow(2, 10));
		ckks.mult(plain_10, temp, g, true); // rescale
		//error
		//if (printmode == "debug")
		//{
		//	res_ctxt = ckks.decode_ctxt(g);
		//	res_ctxt.resize(x0.size());
		//	cout << "Decryption Result of g:" << endl;
		//	printVector(res_ctxt, false);
		//	cout << "Actual calculation of g:" << endl;
		//	printVector(res_g, false);
		//	iteration_scale = int(log2(g.scale()));
		//	if (iteration_scale != ckks.get_scale())
		//	{
		//		cout << "G SCALE ERROR!!!" << endl;
		//		cout << format("iteration scale: {}", iteration_scale) << endl;
		//		cout << format("original scale: {}", ckks.get_scale()) << endl;
		//		return;
		//	}
		//	cout << "--------" << endl;
		//}

		//calculate y = y * h
		ckks.mult(y, h, y, false);
		plain_10 = ckks.encode(1.0, y, pow(2.0, 35));
		ckks.mult(plain_10, y, true); // rescale
		//error
		if (printmode == "debug")
		{
			res_ctxt = ckks.decode_ctxt(y);
			res_ctxt.resize(x0.size());
			res_err = subVectors(res_ctxt, res_y);
			res_err = absVectors(res_err);
			res_real_err = subScalar(res_ctxt, 1/sqrt(x));
			res_real_err = absVectors(res_real_err);

			//cout << "Decryption Result of y:" << endl;
			//printVector(res_ctxt, false);
			//cout << "Actual calculation of y:" << endl;
			//printVector(res_y, false);
			//cout << "Error:" << endl;
			//printVector_10eform(res_err, false);
			//cout << "Error(with real invsqrt):" << endl;
			//printVector_10eform(res_real_err, false);
			cout << "Approximate success?(Plaintext).. \t";
			vector<double> res_err_plain = subScalar(res_y, 1 / sqrt(x));
			res_err_plain = absVectors(res_err_plain);
			printVector_OX(res_err_plain, 1e-4);
			cout << "Approximate success?(Ciphertext)\t";
			printVector_OX(res_real_err, 1e-4);
			iteration_scale = int(log2(y.scale()));
			if (iteration_scale != ckks.get_scale())
			{
				cout << "Y SCALE ERROR!!!" << endl;
				cout << format("iteration scale: {}", iteration_scale) << endl;
				cout << format("original scale: {}", ckks.get_scale()) << endl;
				return;
			}
			cout << "--------" << endl;
		}

		cout << format("Remain Level: {}", y.coeff_modulus_size()) << endl;
		cout << "------------------------------" << endl;

		//zzap-bootstrapping
		if (y.coeff_modulus_size() <= 2 && i != iter - 1)
		{
			cout << "Re-encrypting..." << endl;
			y = ckks.encrypt(res_ctxt);
			g = ckks.encrypt(res_g);
			h = ckks.encrypt(res_h);
			cout << "------------------------------" << endl;
		}
	}
}

void beta_seal(double P, double a, double b, Ciphertext x, ckks_build& ckks)
{
	Plaintext minus1 = ckks.encode(-1, x, 1.0);
	Ciphertext x2;
	ckks.mult(minus1, x, x2, false); //25
	Plaintext div = ckks.encode(1/(b-a), x);
	ckks.mult(div, x2, x2, false); // 50
	Plaintext x1 = ckks.encode(P/(b-a), x2);
	Ciphertext x_sgn;
	ckks.add(x1, x2, x_sgn); // 50
	Plaintext plus1 = ckks.encode(1, x_sgn);
	ckks.add(plus1, x_sgn); // 50
	Plaintext half = ckks.encode(0.5, x_sgn);
	ckks.mult(half, x_sgn, false); // 75
	Plaintext one = ckks.encode(1.0, x_sgn, pow(2, 10));
	ckks.mult(one, x_sgn, true);
	
	// 0.5 * (1+sgn(x1-x2))
	
}