#pragma once

#include "SEAL_MM.h"

class ckks_build {
private:
    unique_ptr<EncryptionParameters> parms;
    vector<int> modulus;
    unique_ptr<SEALContext> context;
    unique_ptr<KeyGenerator> keygen;
    SecretKey sk;
    PublicKey pk;
    RelinKeys rlk;
    GaloisKeys glk;

    unique_ptr<Encryptor> enc;
    unique_ptr<Evaluator> eva;
    unique_ptr<Decryptor> dec;
    unique_ptr<CKKSEncoder> encoder;

    int poly_modulus_degree;
    double scale;

public:
    ckks_build(double scale, size_t pmd);

    void createGaloisKeys(int start, int end);

    Plaintext encode(double input);
    Plaintext encode(const double& input, Ciphertext& ctxt);
    Plaintext encode(const double& input, Ciphertext& ctxt, double scale);
    Plaintext encode(const vector<double>& input);
    Plaintext encode(const vector<double>& input, Ciphertext& ctxt);
    Ciphertext encrypt(const Plaintext& plain);
    Ciphertext encrypt(double input);
    Ciphertext encrypt(const double& input, Ciphertext& ctxt);
    Ciphertext encrypt(const vector<double>& input);

    Plaintext decrypt(Ciphertext& ctxt);
    vector<double> decode(const Plaintext& ptxt);
    vector<double> decode_ctxt(Ciphertext& ctxt);

    void modulus_switch(Plaintext& ptxt, const parms_id_type parms_id);
    void modulus_switch(Ciphertext& ctxt, const parms_id_type parms_id);
    void add(Ciphertext& ctxt1, Ciphertext& ctxt2, Ciphertext& result);
    void add(Ciphertext& ctxt1, Ciphertext& ctxt2);
    void add(Plaintext& ptxt, Ciphertext& ctxt);
    void mult(Ciphertext& ctxt1, Ciphertext& ctxt2, Ciphertext& result);
    void mult(Ciphertext& ctxt1, Ciphertext& ctxt2);
    void mult(Plaintext& ptxt, Ciphertext& ctxt);
    void square(Ciphertext& ctxt);

    Ciphertext exp(const Ciphertext& x, int d);
    //void modulus_equal(Ciphertext& ctxt1, Ciphertext& ctxt2);
    void scale_equal(Ciphertext& ctxt1, Ciphertext& ctxt2);
    void scale_equal(Plaintext& ptxt, Ciphertext& ctxt);

    //void mul_plain_double(Plaintext& ptxt, Ciphertext& ctxt, Ciphertext& destination);
    //void mul_cipher_double(Ciphertext& ctxt1, Ciphertext& ctxt2, Ciphertext& ctxt3, Ciphertext& destination);
    //void add_cipher(Ciphertext& ctxt1, Ciphertext& ctxt2, Ciphertext& destination);
    //void temp_d3_doubleScale(vector<double>& poly, Ciphertext& x, Ciphertext& destination);
    //void evaluate_function_tripleScale(vector<double>& poly, Ciphertext& x, Ciphertext& destination);
    void evaluate_function_tripleScale_v2(vector<double>& poly, Ciphertext& x, Ciphertext& destination);


    //#######################################################
    
    vector<Plaintext> encode_matrix(vector<vector<double>> U, Ciphertext& x);
    Ciphertext rotate(Ciphertext& x, int k);
    Ciphertext linear_transform_sub(Ciphertext& x, Plaintext& coeffs, int rot, int xLength);
    Ciphertext linear_transform(Ciphertext& x, vector<vector<double>>& dVecs, int xLength);
    Ciphertext matrix_multiplication(Ciphertext& A, Ciphertext& B, int originalD);

    void debug_printMatrix(Ciphertext& A, int originalD, string title);
    void debug_vector(vector<double> A, int size, string title);
};
