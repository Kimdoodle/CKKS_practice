#include "header/SEAL_VS.h"

//Generator
ckks_build::ckks_build(string mode, int alpha, int n, int d, int big_moduli, int small_moduli, double scale, size_t pmd)
{
    this->mode = mode;
    this->alpha = alpha;
    this->scale = scale;
    
    parms = make_unique<EncryptionParameters>(scheme_type::ckks);
    parms->set_poly_modulus_degree(pmd);
    
    vector<int> modulus;
    modulus.push_back(big_moduli);
    for (int i = 0; i < d * (2 * n + 2); i++)
        modulus.push_back(small_moduli);
    modulus.push_back(big_moduli);

    int r = 0;
    for (int m : modulus) { r += m; }
    if (CoeffModulus::MaxBitCount(pmd) < r) {
        cerr << "Error: " << CoeffModulus::MaxBitCount(pmd) << " < " << r << endl;
        throw invalid_argument("MaxBitCount Error!");
    }

    parms->set_coeff_modulus(CoeffModulus::Create(pmd, modulus));

    context = make_unique<SEALContext>(*parms);
    print_parameters(*context);
    keygen = make_unique<KeyGenerator>(*context);
    sk = keygen->secret_key();
    keygen->create_public_key(pk);
    keygen->create_relin_keys(rlk);

    enc = make_unique<Encryptor>(*context, pk);
    eva = make_unique<Evaluator>(*context);
    dec = make_unique<Decryptor>(*context, sk);
    encoder = make_unique<CKKSEncoder>(*context);

    //debug - calculate scale factors
    if(mode == "debug")
        calscales();
}

//scale factor debug
void ckks_build::calscales() 
{
    Ciphertext x = encrypt(1.0);
    Ciphertext temp = x;
    scales.push_back(x.scale());
    while (temp.coeff_modulus_size() != 1) {
        mult(temp, x);
        scales.push_back(temp.scale());
    }
}

// encode
Plaintext ckks_build::encode(const vector<double>& input)
{
    Plaintext plain;
    encoder->encode(input, scale, plain);
    return plain;
}
Plaintext ckks_build::encode(const double& input, Ciphertext& ctxt)
{
    Plaintext plain;
    encoder->encode({ input }, ctxt.parms_id(), ctxt.scale(), plain);
    return plain;
}

// encrypt
Ciphertext ckks_build::encrypt(const Plaintext& plain)
{
    Ciphertext ctxt;
    enc->encrypt(plain, ctxt);
    return ctxt;
}

// encode + encrypt
Ciphertext ckks_build::encrypt(double input)
{
    return encrypt(encode({ input }));
}

Ciphertext ckks_build::encrypt(const vector<double>& input)
{
    return encrypt(encode(input));
}

Ciphertext ckks_build::encrypt(const double& input, Ciphertext& ctxt)
{
    return encrypt(encode(input, ctxt));
}

// decrypt
Plaintext ckks_build::decrypt(Ciphertext& ctxt)
{
    Plaintext dec_plain;
    dec->decrypt(ctxt, dec_plain);
    return dec_plain;
}

// decode
vector<double> ckks_build::decode(const Plaintext& ptxt)
{
    vector<double> output;
    encoder->decode(ptxt, output);
    return output;
}

// decrypt + decode
vector<double> ckks_build::decode_ctxt(Ciphertext& ctxt)
{
    return decode(decrypt(ctxt));
}

// Modulus Switch
void ckks_build::modulus_switch(Plaintext& ptxt, const parms_id_type parms_id)
{
    eva->mod_switch_to_inplace(ptxt, parms_id);
}

void ckks_build::modulus_switch(Ciphertext& ctxt, const parms_id_type parms_id)
{
    eva->mod_reduce_to_inplace(ctxt, parms_id);
}

// Add
void ckks_build::add(Ciphertext& ctxt1, Ciphertext& ctxt2, Ciphertext& result)
{
    scale_equal(ctxt1, ctxt2);
    scale_equal(ctxt1, result);
    eva->add(ctxt1, ctxt2, result);
}
void ckks_build::add(Ciphertext& ctxt1, Ciphertext& ctxt2)
{
    scale_equal(ctxt1, ctxt2);
    eva->add_inplace(ctxt1, ctxt2);
}
void ckks_build::add(Plaintext& ptxt, Ciphertext& ctxt)
{
    scale_equal(ptxt, ctxt);
    eva->add_plain_inplace(ctxt, ptxt);
}

// Multiply
void ckks_build::mult(Ciphertext& ctxt1, Ciphertext& ctxt2, Ciphertext& result)
{
    modulus_equal(ctxt1, ctxt2);
    eva->multiply(ctxt1, ctxt2, result);
    eva->relinearize_inplace(result, rlk);
    eva->rescale_to_next_inplace(result);
}
void ckks_build::mult(Ciphertext& ctxt1, Ciphertext& ctxt2)
{
    modulus_equal(ctxt1, ctxt2);
    eva->multiply_inplace(ctxt1, ctxt2);
    eva->relinearize_inplace(ctxt1, rlk);
    eva->rescale_to_next_inplace(ctxt1);
}
void ckks_build::mult(Plaintext& ptxt, Ciphertext& ctxt)
{
    modulus_switch(ptxt, ctxt.parms_id());
    eva->multiply_plain_inplace(ctxt, ptxt);
    eva->rescale_to_next_inplace(ctxt);
}
void ckks_build::square(Ciphertext& ctxt)
{
    eva->square_inplace(ctxt);
    eva->relinearize_inplace(ctxt, rlk);
    eva->rescale_to_next_inplace(ctxt);
}
// exp
Ciphertext ckks_build::exp(const Ciphertext& x, int d)
{
    Ciphertext result, squareX;
    if (d % 2 == 0) {
        result = encrypt(1.0);
    }
    else {
        result = x;
    }

    squareX = x;
    d /= 2;
    while (d > 0) {
        square(squareX);
        if (d % 2) {
            mult(result, squareX);
        }
        d /= 2;
    }
    return result;
}

// 두 암호문의 modulus(level) 통일
void ckks_build::modulus_equal(Ciphertext& ctxt1, Ciphertext& ctxt2)
{
    auto level1 = ctxt1.coeff_modulus_size();
    auto level2 = ctxt2.coeff_modulus_size();

    if (level1 < level2) {
        modulus_switch(ctxt2, ctxt1.parms_id());
    }
    else if (level1 > level2) {
        modulus_switch(ctxt1, ctxt2.parms_id());
    }
}

//두 암호문의 scale 통일
void ckks_build::scale_equal(Ciphertext& ctxt1, Ciphertext& ctxt2)
{
    //while (ctxt1.coeff_modulus_size() != ctxt2.coeff_modulus_size()) {
    //    if (ctxt1.coeff_modulus_size() > ctxt2.coeff_modulus_size()) {
    //        Plaintext p = encode({ 1.0 });
    //        mult(p, ctxt1);
    //    }
    //    else if (ctxt1.coeff_modulus_size() < ctxt2.coeff_modulus_size()) {
    //        Plaintext p = encode({ 1.0 });
    //        mult(p, ctxt2);
    //    }
    //    //debug
    //    if (mode == "debug") {
    //        if(ctxt1.scale() != scales[scales.size() - ctxt1.coeff_modulus_size()])
    //            cout << "RESCALING ERROR!!!! ctxt1." << endl;
    //        if (ctxt2.scale() != scales[scales.size() - ctxt2.coeff_modulus_size()])
    //            cout << "RESCALING ERROR!!!! ctxt2." << endl;
    //    }
    //}
    if (ctxt1.coeff_modulus_size() > ctxt2.coeff_modulus_size()) {
        modulus_switch(ctxt1, ctxt2.parms_id());
        eva->rescale_to_inplace(ctxt1, ctxt2.parms_id());
    }
    else if (ctxt1.coeff_modulus_size() < ctxt2.coeff_modulus_size()) {
        modulus_switch(ctxt2, ctxt1.parms_id());
        eva->rescale_to_inplace(ctxt2, ctxt1.parms_id());
    }
}
void ckks_build::scale_equal(Plaintext& ptxt, Ciphertext& ctxt)
{
    //Todo if necessary.
}
