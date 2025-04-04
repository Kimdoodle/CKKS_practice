#include "header/SEAL_VS.h"

// print vector elements. Choose vector is function's coeff or not.
void printVector(vector<double>& coeffs, bool asFunction, int pre) {
    cout << fixed << setprecision(pre);
    size_t size = coeffs.size();
    for (size_t i = 0; i < size; i++) {
        //if (coeffs[i] == 0.0) continue;

        cout << coeffs[i];

        if (asFunction) {
            cout << "(x" << i << ")";
        }
        if (i != (size - 1)) {
            if (asFunction)
                cout << " + ";
            else
                cout << "\t";
        }
    }
    cout << endl;
}

void print_parameters(const SEALContext& context)
{
    auto& context_data = *context.key_context_data();
    string scheme_name = "CKKS";

    cout << "| Encryption parameters :" << endl;
    cout << "|   scheme: " << scheme_name << endl;
    cout << "|   poly_modulus_degree: " << context_data.parms().poly_modulus_degree() << endl;

    /*
    Print the size of the true (product) coefficient modulus.
    */
    cout << "|   coeff_modulus size: ";
    cout << context_data.total_coeff_modulus_bit_count() << " (";
    auto coeff_modulus = context_data.parms().coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();
    for (size_t i = 0; i < coeff_modulus_size - 1; i++)
    {
        cout << coeff_modulus[i].bit_count() << " + ";
    }
    cout << coeff_modulus.back().bit_count();
    cout << ") bits" << endl;
    cout << "---" << endl;
}

/*
    use between evaluation levels.
    print remaining level, decryption result, real calculation result.
*/
void printStep(Ciphertext& y, vector<double>& real_result, vector<double>& poly, int iter, int size, ckks_build& ckks, int pre)
{
    cout << "d = " << iter + 1 << endl;
    //Remaining Levels
    cout << "\tRemaining Level: " << y.coeff_modulus_size() << endl;
    real_result = polypolyEvaluate(poly, real_result);
    //printResult(y, real_result, size, ckks, pre);
}

/*
    print Error.
*/
void printResult(vector<double>& dec_result, vector<double>& real_result, int pre)
{
    cout << "Difference:\n\t";
    vector<double> error(dec_result.size());
    for (int i = 0; i < dec_result.size(); i++)
        error[i] = abs(dec_result[i] - real_result[i]);
    printVector(error, false, 6);
}