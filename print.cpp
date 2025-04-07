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

void printStep(vector<double>& realValue, vector<double>& poly, vector<double>& fnDec, vector<double> raw_inputs, Ciphertext& y, ckks_build& ckks, string function, int i)
{
    realValue = polypolyEvaluate(poly, realValue);
    fnDec = ckks.decode_ctxt(y);
    fnDec.resize(raw_inputs.size());
    cout << function << "^" << i << " 실제 계산값:\n\t";
    printVector(realValue, false, 6);
    cout << function << "^" << i << " 복호화값:\n\t";
    printVector(fnDec, false, 6);
    cout << "Remaining Levels: " << y.coeff_modulus_size() << endl;
    cout << "Current Scale: " << y.scale() << endl;
    cout << "---" << endl;
}