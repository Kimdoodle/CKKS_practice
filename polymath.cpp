#include "header/SEAL_MM.h"

// 팩토리얼
int factorial(int a, int b) {
    int res = 1;
    for(int i=a; i>b; i--)
        res *= i;
    return res;
}

//밑이 2인 로그
double log2(double x, double base) {
    return log(x)/log(base);
}

// 미분
vector<double> differentiate(vector<double> poly)
{
    vector<double> result;
    for (int i = 1; i < poly.size(); i++) {
        result.push_back(poly[i] * i);
    }
    return result;
}

/*
    sample data in [min, -epsilon] U [epsilon, max]
    repeat iter times.
*/
vector<double> sample_data(double min, double max, double epsilon,  int iter)
{
    vector<double> samples;
    samples.reserve(iter);
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist1(-1.0, -epsilon);
    uniform_real_distribution<double> dist2(epsilon, 1.0);
    for (int i = 0; i < iter; i++) {
        double sample = (rand() % 2 == 0) ? dist1(gen) : dist2(gen);
        samples.push_back(sample);
    }
    return samples;
}

//특정 값을 size만큼 복제한 벡터 생성
vector<double> duplicate_vector(double input, int size)
{
    vector<double> res;
    for (int i = 0; i < size; i++) res.push_back(input);
    return res;
}

//벡터의 모든 요소에 plain 곱셈
vector<double> multPlainPolynomial(vector<double>& v, double scalar)
{
    vector<double> result(v.size());
    for (int i = 0; i < v.size(); i++)
        result[i] = scalar * v[i];
    return result;
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

// (x,y)점들로 다항식을 계산하는 함수(라그랑주 다항식)
vector<double> calculatePoly(const vector<double>& x, const vector<double>& y) {
    int n = x.size();
    vector<double> result(n, 0.0);

    for (int i = 0; i < n; i++) {
        double xi = x[i];
        double yi = y[i];

        vector<double> term = { 1.0 }; // L_i(x) = 1 초기화
        double denominator = 1.0;

        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            double xj = x[j];
            vector<double> poly_term = { -xj, 1.0 }; // (x - x_j)
            term = multPolynomial(term, poly_term);
            denominator *= (xi - xj);
        }

        term = multPlainPolynomial(term, yi / denominator); // y_i / L_i(xi)

        // 결과 다항식에 더하기
        for (int k = 0; k < term.size(); ++k) {
            result[k] += term[k];
        }
    }
    return result;
}

//행렬 전체 프로세스
vector<vector<double>> make_matrix(int inputsize, const string& type)
{
    auto fresh_matrix = type_matrix(inputsize, type);
    //auto p_r_d_matrix = pad_matrix(fresh_matrix, outputsize);
    auto d_matrix = diagonal_matrix(fresh_matrix);
    return d_matrix;
}

// 타입에 맞는 행렬 생성
vector<vector<double>> type_matrix(int originSize, const string& type)
{
    int resultSize = originSize*originSize;
    vector<vector<double>> result(resultSize, vector<double>(resultSize, 0));

    int row, col = 0;
    for (int i = 0; i < originSize; i++)
    {
        for (int j = 0; j < originSize; j++)
        {
            row = originSize * i + j;
            if (type == "sigma")
            {
                col = originSize * i + (i + j) % originSize;
            }
            else if (type == "tau")
            {
                col = originSize * ((i + j) % originSize) + j;
            }
            else
            {
                int k = stoi(type.substr(3));
                string type2 = type.substr(0, 3);
                if (type2 == "phi")
                {
                    col = originSize * i + ((j + k) % originSize);
                }
                else if (type2 == "psi")
                {
                    col = originSize * ((i + k) % originSize) + j;
                }
            }
            result[row][col] = 1;
        }
    }
    return result;
}

// 행렬 U를 dxd행렬로 확장(0.0으로 채움)
vector<vector<double>> pad_matrix(vector<vector<double>> U, int d)
{
    vector<vector<double>> pMatrix;
    int rows = U.size();

    // 실제 데이터 행 처리
    for (const auto& row : U) {
        vector<double> padded_row = row;
        padded_row.resize(d, 0.0);
        pMatrix.push_back(padded_row);
    }

    // 부족한 행을 0으로 padding
    vector<double> zero_row(d, 0.0);
    for (int i = rows; i < d; ++i) {
        pMatrix.push_back(zero_row);
    }

    return pMatrix;
}

// 행렬 U의 padding을 제거 -> dxd 행렬로 변형
vector<vector<double>> ipad_matrix(vector<vector<double>> U, int d)
{
    vector<vector<double>> nopad_U;
    for (int i = 0; i < d; i++) {
        vector<double> row;
        for (int j = 0; j < d; j++) {
            row.push_back(U[i][j]);
        }
        nopad_U.push_back(row);
    }
    return nopad_U;
}

//행렬 U를 (d^2)x1 벡터로 평탄화
vector<double> flatten_matrix(vector<vector<double>> U)
{
    vector<double> flatMatrix;
    for(vector<double> row: U)
        flatMatrix.insert(flatMatrix.end(), row.begin(), row.end());
    return flatMatrix;
}

//벡터를 dxd행렬로 다시 변형
vector<vector<double>> unflatten_matrix(vector<double> U, int colsize)
{
    int rowsize = U.size() / colsize;
    vector<vector<double>> matrix(rowsize, vector<double>(colsize, 0.0));

    for (int i = 0; i < rowsize; i++) {
        for (int j = 0; j < colsize; j++) {
            matrix[i][j] = U[colsize * i + j];
        }
    }
    return matrix;
}

// 행렬에 따른 대각벡터 생성
vector<vector<double>> diagonal_matrix(vector<vector<double>> U)
{
    size_t n = U.size();
    vector<vector<double>> dVecs;

    for (size_t i = 0; i < n; ++i)
    {
        vector<double> dVec;
        for (size_t j = 0; j < n; ++j)
        {
            dVec.push_back(U[j][(j + i) % n]);
        }
        dVecs.push_back(dVec);
    }

    return dVecs;
}

map<int, vector<double>> diagonal_vector(vector<vector<double>> U, string type)
{
    map<int, vector<double>> result;
    size_t len = U.size() * U[0].size();
    int d = static_cast<int>(sqrt(len));

    // Flatten
    vector<double> ct;
    for (const auto& row : U)
        ct.insert(ct.end(), row.begin(), row.end());

    if (type == "sigma")
    {
        for (int k = -d + 1; k < 0; ++k)
        {
            vector<double> uk(len, 0);
            for (int l = 0; l < len; ++l)
            {
                int cond = l - (d + k) * d;
                if (-k <= cond && cond < d)
                    uk[l] = 1;
            }
            result[k] = uk;
        }
        for (int k = 0; k < d; ++k)
        {
            vector<double> uk(len, 0);
            for (int l = 0; l < len; ++l)
            {
                int cond = l - d * k;
                if (0 <= cond && cond < (d - k))
                    uk[l] = 1;
            }
            result[k] = uk;
        }
    }
    else if (type == "tau")
    {
        for (int k = 0; k < d; ++k)
        {
            vector<double> uk(len, 0);
            int l = 0;
            while (k + d * l < static_cast<int>(len))
            {
                uk[k + d * l] = 1;
                ++l;
            }
            result[d * k] = uk;
        }
    }
    else if (type == "phi")
    {
        for (int variant = 0; variant < 2; ++variant)
        {
            vector<double> uk(len, 0);
            for (int l = 0; l < static_cast<int>(len); ++l)
            {
                int mod = l % d;
                if (variant == 0 && 0 <= mod && mod < (d))
                    uk[l] = (mod < d) ? 1 : 0;
                else if (variant == 1 && (d - 0) <= mod && mod < d)
                    uk[l] = 1;
            }
            result[(variant == 0) ? 0 : -d] = uk;
        }
    }
    else if (type == "psi")
    {
        vector<double> uk(len, 1);
        result[0] = uk;
    }

    return result;
}