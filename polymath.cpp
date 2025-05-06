#include "header/SEAL_VS.h"

// factorial
int factorial(int a, int b) {
    int res = 1;
    for(int i=a; i>b; i--)
        res *= i;
    return res;
}

//log(base=2)
double log2(double x, double base) {
    return log(x)/log(base);
}

//differentiate coeff vector
vector<double> differentiate(vector<double> poly)
{
    vector<double> result;
    for (int i = 1; i < poly.size(); i++) {
        result.push_back(poly[i] * i);
    }
    return result;
}

/*
    sample data in [min, max]
    repeat iter times.
*/
vector<double> sample_data(double min, double max, int iter)
{
    vector<double> samples;
    samples.reserve(iter);
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(min, max);
    for (int i = 0; i < iter; i++) {
        double sample = dist(gen);
        samples.push_back(sample);
    }
    return samples;
}

//duplicate vector, 
vector<double> duplicate_vector(double input, int size)
{
    vector<double> res;
    for (int i = 0; i < size; i++) res.push_back(input);
    return res;
}

//multiply scalar in each vector
vector<double> multPlainPolynomial(vector<double>& v, double scalar)
{
    vector<double> result(v.size());
    for (int i = 0; i < v.size(); i++)
        result[i] = scalar * v[i];
    return result;
}
vector<double> multScalar(vector<double>& v, double scalar)
{
    vector<double> result(v.size());
    for (int i = 0; i < v.size(); i++)
        result[i] = scalar * v[i];
    return result;
}

//add scalar to vector element-wise.
vector<double> addScalar(vector<double>& a, double scalar)
{
    vector<double> result(a.size());
    for (int i = 0; i < a.size(); i++) {
        result[i] = a[i] + scalar;
    }
    return result;
}

//add two vectors element-wise.
vector<double> addVectors(vector<double>& a, vector<double>& b)
{
    vector<double> result(a.size());
    for (int i = 0; i < a.size(); i++) {
        result[i] = a[i] + b[i];
    }
    return result;
}

//multiply two vectors element-wise.
vector<double> multVectors(vector<double>& a, vector<double>& b)
{
    vector<double> result(a.size());
    for (int i = 0; i < a.size(); i++) {
        result[i] = a[i] * b[i];
    }
    return result;
}

vector<double> subScalar(vector<double>& a, double scalar)
{
    vector<double> result(a.size());
    for (int i = 0; i < a.size(); i++) {
        result[i] = a[i] - scalar;
    }
    return result;
}

vector<double> subVectors(vector<double>& a, vector<double>& b)
{
    vector<double> result(a.size());
    for (int i = 0; i < a.size(); i++) {
        result[i] = a[i] - b[i];
    }
    return result;
}

vector<double> absVectors(vector<double>& a)
{
    vector<double> result(a.size());
    for (int i = 0; i < a.size(); i++) {
        result[i] = abs(a[i]);
    }
    return result;
}

/*
    ?�항??곱셈???�한 계수 벡터 ?�성
    1. Toeplitz ?�렬???�성
    2. ?�렬 * 벡터 결과�?반환
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

// ?�렬�?벡터??�?
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


// ?�항??곱셈 ?�수
vector<double> multPolynomial(const vector<double>& a, const vector<double>& b) {
    int result_size = static_cast<int>(a.size() + b.size() - 1);
    vector<vector<double>> T = createToeplitzMatrix(a, result_size);

    vector<double> extended_b(result_size, 0);
    for (size_t i = 0; i < b.size(); i++) {
        extended_b[i] = b[i];
    }

    return multiplyMatrixVector(T, extended_b);
}

// ?�항??거듭?�곱 ?�수
vector<double> powerPolynomial(const vector<double>& poly, int exponent) {
    vector<double> result = {1};

    for (int i = 0; i < exponent; i++) {
        result = multPolynomial(result, poly);
    }

    return result;
}

//evaluate polynomial.
double polyEvaluate(const vector<double>& poly, double input) 
{
    double result = 0.0;
    for(int i=0; i<poly.size(); i++) {
        result += poly[i] * pow(input, i);
    }
    return result;
}


//evaluate polynomial d times.
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

// (x,y)?�들�??�항?�을 계산?�는 ?�수(?�그?�주 ?�항??
vector<double> calculatePoly(const vector<double>& x, const vector<double>& y) {
    int n = x.size();
    vector<double> result(n, 0.0);

    for (int i = 0; i < n; i++) {
        double xi = x[i];
        double yi = y[i];

        vector<double> term = { 1.0 }; // L_i(x) = 1 초기??
        double denominator = 1.0;

        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            double xj = x[j];
            vector<double> poly_term = { -xj, 1.0 }; // (x - x_j)
            term = multPolynomial(term, poly_term);
            denominator *= (xi - xj);
        }

        term = multPlainPolynomial(term, yi / denominator); // y_i / L_i(xi)

        // 결과 ?�항?�에 ?�하�?
        for (int k = 0; k < term.size(); ++k) {
            result[k] += term[k];
        }
    }
    return result;
}