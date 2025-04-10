#include "header/SEAL_MM.h"

int main()
{
	// param settings
	vector<vector<double>> A = {
		{1, 2, 3},
		{4, 5, 6},
		{7, 8, 9}
	};

	
	int pmd = 32768;
	double scale = pow(2.0, 60);

	ckks_build ckks = ckks_build(scale, pmd);
	vector<vector<double>> newA = pad_matrix(A, int(sqrt(pmd / 2)));
	vector<double> flatA = flatten_matrix(newA);
	Ciphertext ctA = ckks.encrypt(flatA);

	Ciphertext AA = ckks.matrix_multiplication(ctA, ctA, A.size());
	
	vector<double> decA = ckks.decode_ctxt(AA);
	vector<vector<double>> unflat_decA = unflatten_matrix(decA, A.size());
	printMatrix(unflat_decA, 3);


	//auto time1 = cur_time(); //start time

	//// evaluate
	//Ciphertext y = x;
	//cout << "Initial Level: " << y.coeff_modulus_size() << endl;
	//cout << "Initial Scale: " << y.scale() << endl;
	//cout << "------------------------------------" << endl;


	//auto time2 = cur_time(); // end time
	//calculate_time(time1, time2);

	return 0;
}
