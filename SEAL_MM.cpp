#include "header/SEAL_MM.h"

int main()
{
	// param settings
	//vector<vector<double>> A = {
	//	{1, 2, 3},
	//	{4, 5, 6},
	//	{7, 8, 9}
	//};

	vector<double> A = { 1,2,3,4,5,6,7,8,9 };
	int d = A.size();
	int pmd = 32768;
	double scale = pow(2.0, 60);

	ckks_build ckks = ckks_build(scale, pmd);
	Ciphertext ctA = ckks.encrypt(A);

	auto time1 = cur_time(); //start time

	Ciphertext AA = ckks.matrix_multiplication(ctA, ctA, d);
	vector<double> AAA = ckks.decode_ctxt(AA);
	AAA.resize(d);
	cout << "Matrix Multiplication Result:" << endl;
	printVector(AAA, false, 0);

	auto time2 = cur_time(); // end time
	calculate_time(time1, time2);

	return 0;
}
