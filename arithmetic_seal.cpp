//암호문을 다항식으로 평가

#include "header/SEAL_VS.h"

using namespace std;
using namespace seal;

//두 다항식 연산
Ciphertext cal_x1_x2(char mode, Ciphertext x1, Ciphertext x2, const SEALContext& context, Evaluator& eva, RelinKeys& rlk, double scale) 
{
	Ciphertext result;
	//낮은 level을 검색
	auto level1 = x1.coeff_modulus_size();
	auto level2 = x2.coeff_modulus_size();

	if (level1 < level2) {
		eva.mod_reduce_to_inplace(x2, x1.parms_id());
	}
	else if (level1 > level2) {
		eva.mod_reduce_to_inplace(x1, x2.parms_id());
	}

	//검사
	if (x1.parms_id() != x2.parms_id()) {
		cout << "parms_id 다름." << endl;
	}
	if (mode == '+') {
		eva.add(x1, x2, result);
	}
	else if (mode == 'x') {
		eva.multiply(x1, x2, result);
		eva.relinearize_inplace(result, rlk);
		eva.rescale_to_next_inplace(result);
		result.scale() = scale;
	}

	return result;
}

// x의 거듭제곱(홀수차만 사용)
Ciphertext exp_x(Ciphertext x, int d, Evaluator& eva, SEALContext& context, RelinKeys& rlk, double scale)
{
	Ciphertext result = x;
	Ciphertext squareX = x;
	d /= 2;
	while (d > 0) {
		squareX = cal_x1_x2('x', squareX, squareX, context, eva, rlk, scale);
		if (d % 2) {
			result = cal_x1_x2('x', result, squareX, context, eva, rlk, scale);
		}
		d /= 2;
	}
	return result;
}

// x와 상수 곱셈
Ciphertext mult_x_plain(Plaintext coeff, Ciphertext x, double scale, Evaluator &eva)  
{	
	eva.mod_switch_to_inplace(coeff, x.parms_id());
	eva.multiply_plain_inplace(x, coeff);
	eva.rescale_to_next_inplace(x);
	x.scale() = scale;

	return x;
}


