#include "header/SEAL_VS.h"

/*
	Todo
	1. d를 n에 대한 식으로 변경
	2. epsilon변수 추가 -> 샘플링 시 [-e,e]구간 제외
	3. print.cpp 리팩토링
	4. alpha=8로 설정(정확도 파라미터), scale값과 분리
	5. 계산시간 추가
	6. Rescaling -> 2번당 1번으로 줄임?
*/
int main()
{
	int num = 5; // sample수
	int n = 1; // f_n
	int d = 2; // f_n^d
	int alpha = 23;
	int big_moduli = 30;
	int small_moduli = 23;
	double scale = pow(2.0, alpha);
	size_t pmd = 16384;

	string mode = "debug";
	ckks_build ckks = ckks_build(mode, alpha, n, d, big_moduli, small_moduli, scale, pmd);

	//Sample data
	vector<double> raw_inputs;
	raw_inputs.reserve(num);
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<double> dist(-1.0, 1.0);
	for (int i = 0; i < num; i++) {
		double sample = dist(gen);
		raw_inputs.push_back(sample);
	}
	cout << "sampled x\n\t";
	printVector(raw_inputs, false, 6);
	cout << "---" << endl;

	//f_n 계산
	vector<double> poly = computeF(n);
	cout << "f_" << n << ":\t";
	printVector(poly, true);
	cout << "---" << endl;

	//계산할 함수 선택
	int function_mode;
	cout << "1.sgn(x), 2.abs(x), 3.min/max(x)" << endl;
	cin >> function_mode;
	cout << "---" << endl;
	
	vector<double> realValue;
	realValue.resize(raw_inputs.size());
	Ciphertext result;
	Ciphertext x = ckks.encrypt(raw_inputs);

	switch (function_mode) {
	case 1: //1. sgn(x) 계산
		for (int i = 0; i < raw_inputs.size(); i++) {
			realValue[i] = polypolyEvaluate(poly, raw_inputs[i], d);
		}
		result = sgn_seal(mode, x, poly, raw_inputs, d, ckks);
		printResult(result, realValue, raw_inputs.size(), ckks);
		cout << "---" << endl;
		break;

	case 2: //2. abs(x) 계산
		for (int i = 0; i < raw_inputs.size(); i++) {
			realValue[i] = calAbs(raw_inputs[i], n, d);
		}
		result = abs_seal(mode, x, poly, raw_inputs, d, ckks);
		printResult(result, realValue, raw_inputs.size(), ckks);
		cout << "---" << endl;
		break;

	case 3: //3. max(x) 계산
		double maxValue;
		Ciphertext result;
		if (raw_inputs.size() == 2) {
			maxValue = raw_inputs[0] > raw_inputs[1] ? raw_inputs[0] : raw_inputs[1];
			result = max_seal(mode, x, poly, raw_inputs, d, ckks);
		}
		break;
	}

	////복호화
	//vector<double> plain_outputs = ckks.decode_ctxt(result);
	//plain_outputs.resize(realValue.size());
	//cout << "Decryption Completed." << endl << "---" << endl;

	////결과 출력
	//cout << "암호화 근사 결과" << endl;
	//printVector(plain_outputs, false);
	//cout << "비암호화 근사 결과" << endl;
	//printVector(realValue, false);
	//cout << "오차" << endl;
	//for (int i = 0; i < plain_outputs.size(); i++) {
	//	cout << abs(plain_outputs[i] - realValue[i]) << " , ";
	//}
	//cout << endl;

	

	return 0;
}
