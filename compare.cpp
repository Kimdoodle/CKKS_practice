#include "header/SEAL_VS.h"

using namespace std;


//c_n 계산
double calC(int n) {
    return (2*n+1)/pow(4.0, n)*(factorial(2*n)/(factorial(n)*factorial(n)));
}

//f_n 계산
vector<double> computeF(int n) {
    vector<double> coeff;
    // (2n+2)개의 항
    for(int i=0; i<(2*n+2); i++)
        coeff.push_back(0.0);

    for(int i=0; i<=n; i++) {
        double scalar = 1/pow(4.0, i) * (factorial(2*i)/(factorial(i)*factorial(i)));
        //cout << "scalar= " << scalar << endl;
        vector<double> x = {0, 1};
        vector<double> x2 = {1, 0, -1};
        vector<double> c = multPolynomial(x, powerPolynomial(x2, i));
        //cout << "poly: "; printVector(c, true);
        for(int j=0; j<c.size(); j++)
            coeff[j] += c[j]*scalar;
        //cout << "scalar*poly: "; printVector(coeff, true);
    }
    return coeff;
}

//h_n 계산
vector<double> computeH(int n) {
    vector<double> coeff;
    // 최대차수는 2n+1 -> (2n+2)개의 항
    for(int i=0; i<(2*n+2); i++)
        coeff.push_back(0.0);

    for(int i=0; i<=n; i++) {
        double scalar = factorial(2*i)/(factorial(i)*factorial(i));
        vector<double> x = {-1, 2};
        vector<double> x2 = {0, 1, -1};
        vector<double> c = multPolynomial(x, powerPolynomial(x2, i));
        for(int j=0; j<c.size(); j++)
            coeff[j] += c[j]*scalar;
    }
    return coeff;
}

//sgn(x)
double signFunction(double a, int d) {
    double x = a;
    for(int i=1; i<=d; i++) {
        //vector<double> poly = computeF(i);
        vector<double> poly = computeH(i);
        x = polyEvaluate(poly, (x+1)/2);
    }
    return x;
}

//f_n만을 사용한 비교함수
double newComp(double a, double b, int n, int d) {
    //a, b를 [0,1] 내 데이터로 변환
    while(a>1 || b>1) {
        a /= 2;
        b /= 2;
    }

    double result;
    result = signFunction(a-b, d);
    return (result + 1)/2;
}

//f_n, g_n을 사용한 비교함수
//double newCompG;

////절댓값 함수
//double calAbs(double a, int n, int d) {
//    return a * signFunction(a, d);
//}
//
////최솟값 함수
//double calMin(double a, double b, int n, int d) {
//    return (a+b)/2 - calAbs((a-b), n, d)/2;
//}
//
////최댓값 함수
//double calMax(double a, double b, int n, int d) {
//    return (a+b)/2 + calAbs((a-b), n, d)/2;
//}