
#include "SEAL_VS.h"
#include "arithmetic_seal.h"
#include "print.h"
#include "polymath.h"

using namespace seal;

// sign function
Ciphertext sgn(vector<double> input, vector<double> poly, int scale, int d, string mode,
				CKKSEncoder& encoder, Encryptor& enc, Evaluator& eva, 
				SEALContext& context, RelinKeys& rlk, Decryptor& dec) 
{
	// 1. 암호화
	Plaintext plain_inputs;
	encoder.encode(input, scale, plain_inputs);
	Ciphertext x;
	enc.encrypt(plain_inputs, x);
	cout << "---" << endl;

	// 2. 다항식 평가
	Ciphertext new_x = x;
	Ciphertext temp_y;
	Ciphertext term;
	Ciphertext squaresX;
	Plaintext constant;
	Plaintext coeff;
	vector<double> realValue = input;

	//for debug
	Plaintext temp1;
	vector<double> temp2;

	//Repeat f_n d times
	for (int k = 0; k < d; k++) {
		encoder.encode(poly[1], scale, constant);
		temp_y = mult_x_plain(constant, new_x, scale, eva); // ax^1 계산

		//check ax^1
		if (mode == "debug") {
			cout << "d=" << k + 1 << endl;
			cout << "\t" << poly[1] << "(x)" << endl;
			cout << "\tptxt: ";
			for (int i = 0; i < realValue.size(); i++) {
				cout << realValue[i] * poly[1] << " , ";
			}
			cout << endl;
			dec.decrypt(temp_y, temp1);
			encoder.decode(temp1, temp2);
			temp2.resize(input.size());
			cout << "\tctxt: ";
			printVector(temp2, false);
		}

		//Evaluate polynomial f_n(x)
		for (int i = 3; i < poly.size(); i += 2) {
			encoder.encode(poly[i], new_x.parms_id(), scale, coeff);
			squaresX = exp_x(new_x, i, eva, context, rlk, scale); // x^i
			term = mult_x_plain(coeff, squaresX, scale, eva); // ax^i

			temp_y = cal_x1_x2('+', temp_y, term, context, eva, rlk, scale);

			//check ax^3,5,7...
			if (mode == "debug") {
				cout << "\t--" << endl;
				cout << "\t" << poly[i] << "(x^" << i << "): " << endl;
				cout << "\tptxt: ";
				for (int j = 0; j < realValue.size(); j++) {
					cout << pow(realValue[j], i) * poly[i] << " , ";
				}
				cout << endl;
				dec.decrypt(term, temp1);
				encoder.decode(temp1, temp2);
				temp2.resize(input.size());
				cout << "\tctxt";
				printVector(temp2, false);
			}
		}
		new_x = temp_y;

		//check
		if (mode == "debug") {
			for (int i = 0; i < realValue.size(); i++) {
				realValue[i] = polyEvaluate(poly, realValue[i]);
			}
			Plaintext output_plain;
			dec.decrypt(new_x, output_plain);
			vector<double> output_vector;
			encoder.decode(output_plain, output_vector);
			output_vector.resize(input.size());
			cout << "\t--" << endl;
			cout << "\t실제 값: "; printVector(realValue, false);
			cout << "\t복호화값: "; printVector(output_vector, false);
			cout << "--" << endl;
		}
	}
	return new_x;
}

// abs function
// function calculates all abs value in input.
Ciphertext abs_seal(vector<double> input, vector<double> poly, int scale, int d,
					CKKSEncoder& encoder, Encryptor& enc, Evaluator& eva,
					SEALContext& context, RelinKeys& rlk, Decryptor& dec)
{
	Ciphertext x2;
	Plaintext x2_plain;
	encoder.encode(input, scale, x2_plain);
	enc.encrypt(x2_plain, x2);

	Ciphertext sgnx = sgn(input, poly, scale, d, "debug", encoder, enc, eva, context, rlk, dec);

	return cal_x1_x2('x', x2, sgnx, context, eva, rlk, scale);
}

// min/max function
// input size must be 2.
Ciphertext minMax_seal(vector<double> input, vector<double> poly, int scale, int d, string mode,
						CKKSEncoder& encoder, Encryptor& enc, Evaluator& eva,
						SEALContext& context, RelinKeys& rlk, Decryptor& dec)
{
	Ciphertext result;
	if (input.size() != 2) {
		cout << "Input size error!!" << endl;
		return result;
	}
	// |a+b|/2
	vector<double> ssum = { input[0] + input[1] };
	Ciphertext ssum_abs = abs_seal(ssum, poly, scale, d, encoder, enc, eva, context, rlk, dec);
	Plaintext half;
	encoder.encode({ 0.5 }, scale, half);
	ssum_abs = mult_x_plain(half, ssum_abs, scale, eva);
	
	// |a-b|/2
	vector<double> mmin = { input[0] - input[1] };
	Ciphertext mmin_abs = abs_seal(mmin, poly, scale, d, encoder, enc, eva, context, rlk, dec);
	if(mode == "max")
		encoder.encode({ 0.5 }, scale, half);
	else if(mode == "min")
		encoder.encode({ -0.5 }, scale, half);
	mmin_abs = mult_x_plain(half, mmin_abs, scale, eva);

	// calculate
	result = cal_x1_x2('+', ssum_abs, mmin_abs, context, eva, rlk, scale);
	return result;
}