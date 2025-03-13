#include "SEAL_VS.h"
#include "compare.h"
#include "print.h"
#include "polymath.h"
#include "arithmetic_seal.h"
#include "function_seal.h"

using namespace std;
using namespace seal;

int main()
{
	int num = 5; // sample수
	int n = 1; // f_n
	int d = 4; // f_n^d
	int alpha = 23;
	int big_moduli = 30;
	int small_moduli = 23;
	double scale = pow(2.0, alpha);
	size_t pmd = 16384;

	EncryptionParameters parms(scheme_type::ckks);
	parms.set_poly_modulus_degree(pmd);

	vector<int> modulus;
	modulus.push_back(big_moduli);
	for (int i = 0; i < d*(2 * n + 2); i++)
		modulus.push_back(small_moduli);
	modulus.push_back(big_moduli);
	
	//debug - MaxBitCount
	int r = 0;
	for (int m : modulus) { r += m;	}
	if (CoeffModulus::MaxBitCount(pmd) < r) {
		cout << "MaxBitCount Error!" << endl;
		cout << CoeffModulus::MaxBitCount(pmd) << " < " << r << endl;
		return 0;
	}

	parms.set_coeff_modulus(CoeffModulus::Create(pmd, modulus));

	SEALContext context(parms);
	print_parameters(context);
	KeyGenerator keygen(context);
	auto sk = keygen.secret_key();
	PublicKey pk;
	keygen.create_public_key(pk);
	RelinKeys rlk;
	keygen.create_relin_keys(rlk);

	Encryptor enc(context, pk);
	Evaluator eva(context);
	Decryptor dec(context, sk);
	CKKSEncoder encoder(context);
	size_t slot_count = encoder.slot_count();

	//입력 - 임의로 num개 데이터 샘플링
	vector<double> raw_inputs;
	raw_inputs.reserve(num);
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<double> dist(-1.0, 1.0);
	cout << "샘플링한 x: ";
	for (int i = 0; i < num; i++) {
		double sample = dist(gen);
		raw_inputs.push_back(sample);
		cout << sample << " ";
	}
	cout << endl;
	cout << "---" << endl;

	//f_n 계산
	vector<double> poly = computeF(n);
	cout << "f_n 계산결과: " << endl;
	printVector(poly, true);
	cout << "---" << endl;

	//1. sgn(x) 계산
	vector<double> realValue;
	realValue.reserve(raw_inputs.size());
	for (int i = 0; i < raw_inputs.size(); i++) {
		realValue[i] = polypolyEvaluate(poly, raw_inputs[i], d);
	}
	Ciphertext sgnx = sgn(raw_inputs, poly, scale, d, "debug", encoder, enc, eva, context, rlk, dec);
	cout << "---" << endl;

	//2. abs(x) 계산
	for (int i = 0; i < raw_inputs.size(); i++) {
		realValue[i] = abs(raw_inputs[i]);
	}
	Ciphertext absx = abs_seal(raw_inputs, poly, scale, d, encoder, enc, eva, context, rlk, dec);
	
	//3. min/max(x) 계산
	int minValue, maxValue;
	Ciphertext min_ctxt, max_ctxt;
	if (raw_inputs.size() == 2) {
		int a = raw_inputs[0];
		int b = raw_inputs[1];
		if (a < b) {
			minValue = a;
			maxValue = b;
		}
		else if (a > b) {
			minValue = b;
			maxValue = a;
		}
		min_ctxt = minMax_seal(raw_inputs, poly, scale, d, "min", encoder, enc, eva, context, rlk, dec);
		max_ctxt = minMax_seal(raw_inputs, poly, scale, d, "max", encoder, enc, eva, context, rlk, dec);
	}

	//복호화
	Plaintext plain_outputs;
	cout << "Remaining Level: " << sgnx.coeff_modulus_size() << endl;
	dec.decrypt(sgnx, plain_outputs);
	vector<double> raw_outputs;
	encoder.decode(plain_outputs, raw_outputs);
	raw_outputs.resize(raw_inputs.size());
	cout << "Decryption Completed." << endl << "---" << endl;

	//결과 출력
	cout << "암호화 근사 결과" << endl;
	printVector(raw_outputs, false);
	cout << "비암호화 근사 결과" << endl;
	printVector(realValue, false);
	cout << "오차" << endl;
	for (int i = 0; i < raw_outputs.size(); i++) {
		cout << abs(raw_outputs[i] - realValue[i]) << " , ";
	}
	cout << endl;

	return 0;
}
