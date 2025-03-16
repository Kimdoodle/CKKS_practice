#include "header/SEAL_VS.h"

// 팩토리얼
int factorial(int a) {
    int res = 1;
    for(int i=a; i>0; i--)
        res *= i;
    return res;
}

//밑이 2인 로그
double log2(double x, double base) {
    return log(x)/log(base);
}

//특정 값을 size만큼 복제한 벡터 생성
vector<double> duplicate_vector(double input, int size)
{
    vector<double> res;
    for (int i = 0; i < size; i++) res.push_back(input);
    return res;
}

/*
    다항식 곱셈을 통한 계수 벡터 생성
    1. Toeplitz 행렬을 생성
    2. 행렬 * 벡터 결과를 반환
*/
vector<vector<double>> createToeplitzMatrix(const vector<double>& coeffs, int result_size) {
    vector<vector<double>> T(result_size, vector<double>(result_size, 0));

    size_t coeff_size = coeffs.size();
    for (size_t i = 0; i < coeff_size; i++) {
        for (size_t j = 0; j < result_size - i; j++) {
            T[i + j][j] = coeffs[i];
        }
    }

    return T;
}

// 행렬과 벡터의 곱
vector<double> multiplyMatrixVector(const vector<vector<double>>& matrix, const vector<double>& vec) {
    size_t size = matrix.size();
    vector<double> result(size, 0);

    for (size_t i = 0; i < size; i++) {
        for (size_t j = 0; j < vec.size(); j++) {
            result[i] += matrix[i][j] * vec[j];
        }
    }

    return result;
}


// 다항식 곱셈 함수
vector<double> multPolynomial(const vector<double>& a, const vector<double>& b) {
    int result_size = static_cast<int>(a.size() + b.size() - 1);
    vector<vector<double>> T = createToeplitzMatrix(a, result_size);

    vector<double> extended_b(result_size, 0);
    for (size_t i = 0; i < b.size(); i++) {
        extended_b[i] = b[i];
    }

    return multiplyMatrixVector(T, extended_b);
}

// 다항식 거듭제곱 함수
vector<double> powerPolynomial(const vector<double>& poly, int exponent) {
    vector<double> result = {1};

    for (int i = 0; i < exponent; i++) {
        result = multPolynomial(result, poly);
    }

    return result;
}

//다항식 계산함수
double polyEvaluate(const vector<double>& poly, double input) {
    double result = 0.0;
    for(int i=0; i<poly.size(); i++) {
        result += poly[i] * pow(input, i);
    }
    return result;
}


//다항식 반복 계산함수
double polypolyEvaluate(const vector<double>& poly, double input, int d) {
    double x = input;
    for (int i = 0; i < d; i++) {
        x = polyEvaluate(poly, x);
    }
    return x;
}
vector<double> polypolyEvaluate(const vector<double>& poly, vector<double>& input) {
    vector<double> result;
    result.resize(input.size());
    for (int i = 0; i < input.size(); i++) {
        result[i] = polyEvaluate(poly, input[i]);
    }
    return result;
}