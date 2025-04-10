#include "header/SEAL_MM.h"

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

// g_n 계산 (Remez Algorithm)
vector<double> computeG(int n, double tau, double pre, double a, double b) 
{
    // 1. [a,b] 구간을 동일한 간격으로 나누는 x (n+1)개 설정
    vector<double> allX(n + 1);
    double step = (b - a) / n;
    for (int i = 0; i <= n; i++) {
        allX[i] = a + step * i;
    }

    // 2. p(x)-f(x)가 허용오차 E 이하인 p(x)를 검색 (f(x): y=1)
    vector<double> allY(n + 1);
    for (int i = 0; i <= n; i++) {
        allY[i] = 1 + tau * pow(-1.0, i);
    }

    // 다항식 근사 계산 (p(x))
    vector<double> p = calculatePoly(allX, allY);

    // 3. 새로운 극대점 찾기
    vector<double> newX(n + 2);
    for (int i = 0; i <= n + 1; i++) {
        newX[i] = (allX[i] + allX[i + 1]) / 2.0;
    }

    // 4. 새로운 극대점에서 오차를 다시 계산
    vector<double> error(n + 2);
    for (int i = 0; i <= n + 1; i++) {
        error[i] = abs(polyEvaluate(p, newX[i]) - 1);
    }

    // 5. 최대 오차 위치 찾기
    int maxErrorIndex = 0;
    for (int i = 1; i <= n + 1; i++) {
        if (error[i] > error[maxErrorIndex]) {
            maxErrorIndex = i;
        }
    }

    // 6. 새로운 극대점과 오차를 기준으로 g_n 갱신
    vector<double> g(n + 1);
    for (int i = 0; i <= n; i++) {
        g[i] = allX[i] + pre * (newX[maxErrorIndex] - allX[i]);
    }

    return g;
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

//절댓값 함수
double calAbs(double a, int n, int d) {
    return a * signFunction(a, d);
}

//최솟값 함수
double calMin(double a, double b, int n, int d) {
    return (a+b)/2 - calAbs((a-b), n, d)/2;
}

//최댓값 함수
double calMax(double a, double b, int n, int d) {
    return (a+b)/2 + calAbs((a-b), n, d)/2;
}