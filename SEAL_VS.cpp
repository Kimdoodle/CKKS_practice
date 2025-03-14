#include "header/SEAL_VS.h"

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

	ckks_build ckks = ckks_build("nodebug", alpha, n, d, big_moduli, small_moduli, scale, pmd);

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

	//계산할 함수 선택
	int mode;
	cout << "1.sgn(x), 2.abs(x), 3.min/max(x)" << endl;
	cin >> mode;
	cout << "---" << endl;
	
	vector<double> realValue;
	realValue.resize(raw_inputs.size());
	Ciphertext result;
	Ciphertext x = ckks.encrypt(raw_inputs);
	switch (mode) {
	case 1:
		//1. sgn(x) 계산
		cout << "sgn(x)" << endl;
		for (int i = 0; i < raw_inputs.size(); i++) {
			realValue[i] = polypolyEvaluate(poly, raw_inputs[i], d);
		}
		result = sgn("debug", x, poly, d, ckks);
		cout << "---" << endl;
		break;
	case 2:
		//2. abs(x) 계산
		cout << "abs(x)" << endl;
		for (int i = 0; i < raw_inputs.size(); i++) {
			realValue[i] = calAbs(raw_inputs[i], n, d);
		}
		result = abs_seal(x, poly, d, ckks);
		cout << "---" << endl;
		break;
	case 3:
		//3. max(x) 계산
		double maxValue;
		Ciphertext max_ctxt;
		if (raw_inputs.size() == 2) {
			maxValue = raw_inputs[0] > raw_inputs[1] ? raw_inputs[0] : raw_inputs[1];
			max_ctxt = max_seal(raw_inputs, poly, d, ckks);
		}
		break;
	}

	//복호화
	cout << "Remaining Level: " << result.coeff_modulus_size() << endl;
	vector<double> plain_outputs = ckks.decode_ctxt(result);
	cout << "Decryption Completed." << endl << "---" << endl;

	//결과 출력
	cout << "암호화 근사 결과" << endl;
	printVector(plain_outputs, false);
	cout << "비암호화 근사 결과" << endl;
	printVector(realValue, false);
	cout << "오차" << endl;
	for (int i = 0; i < plain_outputs.size(); i++) {
		cout << abs(plain_outputs[i] - realValue[i]) << " , ";
	}
	cout << endl;

	

	return 0;
}
