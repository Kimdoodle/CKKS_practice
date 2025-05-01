#include "header/SEAL_VS.h"

////calculate c_n
//double calC(int n) {
//    return (2*n+1)/pow(4.0, n)*(factorial(2*n)/(factorial(n)*factorial(n)));
//}
//
////calculate f_n
//vector<double> computeF(int n) {
//    vector<double> coeff;
//    // (2n+2)
//    for(int i=0; i<(2*n+2); i++)
//        coeff.push_back(0.0);
//
//    for(int i=0; i<=n; i++) {
//        double scalar = 1/pow(4.0, i) * (factorial(2*i)/(factorial(i)*factorial(i)));
//        //cout << "scalar= " << scalar << endl;
//        vector<double> x = {0, 1};
//        vector<double> x2 = {1, 0, -1};
//        vector<double> c = multPolynomial(x, powerPolynomial(x2, i));
//        //cout << "poly: "; printVector(c, true);
//        for(int j=0; j<c.size(); j++)
//            coeff[j] += c[j]*scalar;
//        //cout << "scalar*poly: "; printVector(coeff, true);
//    }
//    return coeff;
//}
//
////h_n Í≥ÑÏÇ∞
//vector<double> computeH(int n) {
//    vector<double> coeff;
//    // ÏµúÎ?Ï∞®Ïàò??2n+1 -> (2n+2)Í∞úÏùò ??
//    for(int i=0; i<(2*n+2); i++)
//        coeff.push_back(0.0);
//
//    for(int i=0; i<=n; i++) {
//        double scalar = factorial(2*i)/(factorial(i)*factorial(i));
//        vector<double> x = {-1, 2};
//        vector<double> x2 = {0, 1, -1};
//        vector<double> c = multPolynomial(x, powerPolynomial(x2, i));
//        for(int j=0; j<c.size(); j++)
//            coeff[j] += c[j]*scalar;
//    }
//    return coeff;
//}
//
//// g_n Í≥ÑÏÇ∞ (Remez Algorithm)
//vector<double> computeG(int n, double tau, double pre, double a, double b) 
//{
//    // 1. [a,b] Íµ¨Í∞Ñ???ôÏùº??Í∞ÑÍ≤©?ºÎ°ú ?òÎàÑ??x (n+1)Í∞??§Ï†ï
//    vector<double> allX(n + 1);
//    double step = (b - a) / n;
//    for (int i = 0; i <= n; i++) {
//        allX[i] = a + step * i;
//    }
//
//    // 2. p(x)-f(x)Í∞Ä ?àÏö©?§Ï∞® E ?¥Ìïò??p(x)Î•?Í≤Ä??(f(x): y=1)
//    vector<double> allY(n + 1);
//    for (int i = 0; i <= n; i++) {
//        allY[i] = 1 + tau * pow(-1.0, i);
//    }
//
//    // ?§Ìï≠??Í∑ºÏÇ¨ Í≥ÑÏÇ∞ (p(x))
//    vector<double> p = calculatePoly(allX, allY);
//
//    // 3. ?àÎ°ú??Í∑πÎ???Ï∞æÍ∏∞
//    vector<double> newX(n + 2);
//    for (int i = 0; i <= n + 1; i++) {
//        newX[i] = (allX[i] + allX[i + 1]) / 2.0;
//    }
//
//    // 4. ?àÎ°ú??Í∑πÎ??êÏóê???§Ï∞®Î•??§Ïãú Í≥ÑÏÇ∞
//    vector<double> error(n + 2);
//    for (int i = 0; i <= n + 1; i++) {
//        error[i] = abs(polyEvaluate(p, newX[i]) - 1);
//    }
//
//    // 5. ÏµúÎ? ?§Ï∞® ?ÑÏπò Ï∞æÍ∏∞
//    int maxErrorIndex = 0;
//    for (int i = 1; i <= n + 1; i++) {
//        if (error[i] > error[maxErrorIndex]) {
//            maxErrorIndex = i;
//        }
//    }
//
//    // 6. ?àÎ°ú??Í∑πÎ??êÍ≥º ?§Ï∞®Î•?Í∏∞Ï??ºÎ°ú g_n Í∞±Ïã†
//    vector<double> g(n + 1);
//    for (int i = 0; i <= n; i++) {
//        g[i] = allX[i] + pre * (newX[maxErrorIndex] - allX[i]);
//    }
//
//    return g;
//}
//
////sgn(x)
//double signFunction(double a, int d) {
//    double x = a;
//    for(int i=1; i<=d; i++) {
//        //vector<double> poly = computeF(i);
//        vector<double> poly = computeH(i);
//        x = polyEvaluate(poly, (x+1)/2);
//    }
//    return x;
//}
//
////f_nÎßåÏùÑ ?¨Ïö©??ÎπÑÍµê?®Ïàò
//double newComp(double a, double b, int n, int d) {
//    //a, bÎ•?[0,1] ???∞Ïù¥?∞Î°ú Î≥Ä??
//    while(a>1 || b>1) {
//        a /= 2;
//        b /= 2;
//    }
//
//    double result;
//    result = signFunction(a-b, d);
//    return (result + 1)/2;
//}
//
////f_n, g_n???¨Ïö©??ÎπÑÍµê?®Ïàò
////double newCompG;
//
////?àÎåìÍ∞??®Ïàò
//double calAbs(double a, int n, int d) {
//    return a * signFunction(a, d);
//}
//
////ÏµúÏÜüÍ∞??®Ïàò
//double calMin(double a, double b, int n, int d) {
//    return (a+b)/2 - calAbs((a-b), n, d)/2;
//}
//
////ÏµúÎåìÍ∞??®Ïàò
//double calMax(double a, double b, int n, int d) {
//    return (a+b)/2 + calAbs((a-b), n, d)/2;
//}

//Newton method.
double newton(double x, double y)
{
    return 0.5 * y * (3 - (x * y * y));
}

double iter_newton(double x, double y, int iter)
{
    double y2 = x;
    for (int i = 0; i < iter; i++)
    {
        y2 = newton(y, y2);
    }
    return y2;
}

double calculate_k1(double low, double high, int iter, double delta, double err, string printmode)
{
    int count = 0;
    debug_print(format("Count {}", count), printmode);
    debug_print(format("Range: {:.6f} ~ {:.6f}", low, high), printmode);

    while ((high - low) >= delta)
    {
        count += 1;
        debug_print(format("Count {}", count), printmode);

        double mid = (high + low) / 2;
        double val = abs(iter_newton(mid, 1.0, iter) - 1);
        debug_print(format("val: {}", val), printmode);

        if (val <= err)
        {
            high = mid;
            debug_print("high -> mid", printmode);
        }
        else
        {
            low = mid + delta;
            debug_print("Increasing low(+delta)", printmode);
        }
        debug_print(format("Low: {:.6f} \t High: {:.6f}", low, high), printmode);
    }

    return low;
}

double calculate_k2(double low, double high, int iter, double delta, double err, string printmode)
{
    int count = 0;
    debug_print(format("Count {}", count), printmode);
    debug_print(format("Range: {:.6f} ~ {:.6f}", low, high), printmode);

    while ((high - low) >= delta)
    {
        count += 1;
        debug_print(format("Count {}", count), printmode);

        double mid = (high + low) / 2;
        double val = abs(iter_newton(mid, 1.0, iter) - 1);

        if (val <= err)
        {
            low = mid;
            debug_print("low -> mid", printmode);
        }
        else
        {
            high = mid - delta;
            debug_print("Decreasing high(-delta)", printmode);
        }

        debug_print(format("Low: {:.6f} \t High: {:.6f}", low, high), printmode);
    }

    return high;
}

pair<double, double> find_bounds(int iter, double delta, double err, const string& mode)
{
    debug_print("Find Lower Bound k1.", mode);
    double low1 = 2 * delta - 1;
    double high1 = 1.0;
    double k1 = calculate_k1(low1, high1, iter, delta, err, mode);

    debug_print(string(30, '-'), mode);
    debug_print("Find Upper Bound k2.", mode);
    double low2 = 1.0;
    double high2 = 2 * sqrt(3.0) - 1;
    double k2 = calculate_k2(low2, high2, iter, delta, err, mode);

    debug_print(string(30, '-'), mode);
    return { k1, k2 };
}

// search x2 using newton-raphson method.
double solve_eq5(double k1, double k2, double x, double x0, int max_iter, double tol)
{
    double x2 = x0;
    for (int i = 0; i < max_iter; ++i)
    {
        // f(x2)
        double f = pow(k2, 2) * pow(x, 3)
            - 6 * pow(k2, 2) * pow(x, 2) * x2
            + 9 * pow(k2, 2) * x * pow(x2, 2)
            - 4 * pow(k1, 2) * pow(x2, 3);

        // f'(x2)
        double df = -6 * pow(k2, 2) * pow(x, 2)
            + 18 * pow(k2, 2) * x * x2
            - 12 * pow(k1, 2) * pow(x2, 2);

        double dx = f / df;
        x2 = x2 - dx;

        if (abs(dx) < tol)
            break;
    }
    return x2;
}

double solve_eq5_for_x(double k1, double k2, double x2_fixed, double x0, int max_iter, double tol) 
{
    double x = x0;
    for (int i = 0; i < max_iter; ++i) {
        double f = pow(k2, 2) * pow(x, 3)
            - 6 * pow(k2, 2) * pow(x, 2) * x2_fixed
            + 9 * pow(k2, 2) * x * pow(x2_fixed, 2)
            - 4 * pow(k1, 2) * pow(x2_fixed, 3);

        double df = 3 * pow(k2, 2) * pow(x, 2)
            - 12 * pow(k2, 2) * x * x2_fixed
            + 9 * pow(k2, 2) * pow(x2_fixed, 2);

        double dx = f / df;
        x -= dx;

        if (abs(dx) < tol) break;
    }
    return x;
}

//comparison functions
double approx_sign(double x, int dg, int df)
{
    vector<double> coeff_f = { 0, 35, 0, -35, 0, 21, 0, -5 }; //f3
    coeff_f = multPlainPolynomial(coeff_f, 1/pow(2, 4));
    vector<double> coeff_g = { 0, 4589, 0, -16577, 0, 25614, 0, -12860 }; //g3
    coeff_g = multPlainPolynomial(coeff_g, 1/pow(2, 10));

    double result_g =(polypolyEvaluate(coeff_g, x, dg));
    double result_f = polypolyEvaluate(coeff_f, result_g, df);

    return result_f;
}

double approx_comp(double a, double b, int dg, int df)
{
    double sign = approx_sign(a - b, dg, df);
    return (sign + 1) / 2;
}