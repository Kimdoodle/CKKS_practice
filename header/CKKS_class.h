#ifndef CKKS_CLASS_H
#define CKKS_CLASS_H

#include "SEAL_VS.h"

class ckks_build {
private:
    unique_ptr<EncryptionParameters> parms;
    unique_ptr<SEALContext> context;
    unique_ptr<KeyGenerator> keygen;
    SecretKey sk;
    PublicKey pk;
    RelinKeys rlk;

    unique_ptr<Encryptor> enc;
    unique_ptr<Evaluator> eva;
    unique_ptr<Decryptor> dec;
    unique_ptr<CKKSEncoder> encoder;

    int alpha;
    double scale;
    string mode;
    vector<double> scales;

public:
    ckks_build(string mode, int alpha, int n, int d, int big_moduli, int small_moduli, double scale, size_t pmd);
    void calscales();

    Plaintext encode(const vector<double>& input);
    Plaintext encode(const double& input, Ciphertext& ctxt);
    Ciphertext encrypt(const Plaintext& plain);
    Ciphertext encrypt(double input);
    Ciphertext encrypt(const vector<double>& input);
    Ciphertext encrypt(const double& input, Ciphertext& ctxt);

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
    void modulus_equal(Ciphertext& ctxt1, Ciphertext& ctxt2);
    void scale_equal(Ciphertext& ctxt1, Ciphertext& ctxt2);
    void scale_equal(Plaintext& ptxt, Ciphertext& ctxt);
};

#endif // CKKS_CLASS_H
