#include "header/SEAL_VS.h"

//Generator
ckks_build::ckks_build(string mode, string scaleMode, int n, int d, int big_moduli, int small_moduli, double scale, size_t pmd)
{
    this->mode = mode;
    this->scale = scale;
    
    parms = make_unique<EncryptionParameters>(scheme_type::ckks);
    parms->set_poly_modulus_degree(pmd);
    

    if (scaleMode == "single") {    //single Scale
        modulus_chain_mode1(big_moduli, small_moduli, d);
    }
    else if (scaleMode == "double") {    //double Scale
        modulus = { 59, 41, 41 };
        modulus_chain_mode3(big_moduli, small_moduli, d);
    }
    else if (scaleMode == "triple") {    //triple Scale
        //modulus = { 60, 60, 60, 47 };
        modulus_chain_mode2(big_moduli, small_moduli, 2, d);
    }

    //check MaxBitCount
    int r = 0;
    for (int m : modulus) { r += m; }
    if (CoeffModulus::MaxBitCount(pmd) < r) {
        cout << "Error: Max Modulus Size " << CoeffModulus::MaxBitCount(pmd) << " < " << r << endl;
        exit(0);
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
    if(mode == "debug_scale")
        calscales();
}

/*  
    make modulus chain. 
    Mode1: big moduli at front/end, small moduli fill rest. 
    Total size: iter+2.
*/
void ckks_build::modulus_chain_mode1(int big_moduli, int small_moduli, int iter) 
{
    modulus.push_back(big_moduli);
    for (int i = 0; i < iter; i++)
        modulus.push_back(small_moduli);
    modulus.push_back(big_moduli);
}

/*
    make modulus chain.
    Mode 2: big moduli placed iter1 times at front, and once at the end. Small moduli fill rest.
    Total size : iter1 + iter2 + 1.
*/
void ckks_build::modulus_chain_mode2(int big_moduli, int small_moduli, int iter1, int iter2)
{
    for(int i=0; i<iter1; i++)
        modulus.push_back(big_moduli);
    for (int i = 0; i < iter2; i++)
        modulus.push_back(small_moduli);
    modulus.push_back(big_moduli);
}

/*
    make modulus chain.
    Mode 3: base modulus already set. push small_moduli iter times, big_moduli at the end.
*/
void ckks_build::modulus_chain_mode3(int big_moduli, int small_moduli, int iter1)
{
    for (int i = 0; i < iter1; i++)
        modulus.push_back(small_moduli);
    modulus.push_back(big_moduli);
}

//scale factor debug
void ckks_build::calscales() 
{
    int d = 1;
    Ciphertext x = encrypt(1.0);
    Ciphertext temp = x;
    scales.push_back(x.scale());
    while (temp.coeff_modulus_size() != 1) {
        mult(temp, x);
        scales.push_back(temp.scale());
        cout << "Level " << d << ":\t" << temp.scale() << endl;
        d++;
    }
}

// encode coeff. Only 1 value needed.
Plaintext ckks_build::encode(double input)
{
    Plaintext plain;
    encoder->encode(input, scale, plain);
    return plain;
}
// encode coeff, but scale, parms_id will be equal to ctxt's.
Plaintext ckks_build::encode(const double& input, Ciphertext& ctxt)
{
    Plaintext plain;
    encoder->encode(input, ctxt.parms_id(), ctxt.scale(), plain);
    return plain;
}
//encode coeff, but only parms_id will be equal to ctxt's.
Plaintext ckks_build::encode(const double& input, Ciphertext& ctxt, double scale)
{
    Plaintext plain;
    encoder->encode(input, ctxt.parms_id(), scale, plain);
    return plain;
}

// encode input value. vector needed.
Plaintext ckks_build::encode(const vector<double>& input)
{
    Plaintext plain;
    encoder->encode(input, scale, plain);
    return plain;
}

// encrypt plaintext.
Ciphertext ckks_build::encrypt(const Plaintext& plain)
{
    Ciphertext ctxt;
    enc->encrypt(plain, ctxt);
    return ctxt;
}

// encode+encrypt plain value.
Ciphertext ckks_build::encrypt(double input)
{
    return encrypt(encode(input));
}

// encode+encrypt plain value, but params will be equal with ctxt's.
Ciphertext ckks_build::encrypt(const double& input, Ciphertext& ctxt)
{
    return encrypt(encode(input, ctxt));
}

// encode+encrypt plain vector
Ciphertext ckks_build::encrypt(const vector<double>& input)
{
    return encrypt(encode(input));
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
    //scale_equal(ptxt, ctxt);
    eva->add_plain_inplace(ctxt, ptxt);
}

// Multiply
void ckks_build::mult(Ciphertext& ctxt1, Ciphertext& ctxt2, Ciphertext& result)
{
    //modulus_equal(ctxt1, ctxt2);
    scale_equal(ctxt1, ctxt2);
    eva->multiply(ctxt1, ctxt2, result);
    eva->relinearize_inplace(result, rlk);
    eva->rescale_to_next_inplace(result);
}

//Multiply 2 ciphertexts.
void ckks_build::mult(Ciphertext& ctxt1, Ciphertext& ctxt2)
{
    scale_equal(ctxt1, ctxt2);
    eva->multiply_inplace(ctxt1, ctxt2);
    eva->relinearize_inplace(ctxt1, rlk);
    eva->rescale_to_next_inplace(ctxt1);
}

//Multiply plaintext / ciphertext. two data's scale, modulus_level should be equal.
void ckks_build::mult(Plaintext& ptxt, Ciphertext& ctxt)
{
    eva->multiply_plain_inplace(ctxt, ptxt);
    eva->rescale_to_next_inplace(ctxt);
}

//square ciphertext.
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
    result = encrypt(1.0);

    squareX = x;
    while(d > 0) {
        if (d % 2) {
            mult(result, squareX);
        }
        square(squareX);
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
    while (ctxt1.coeff_modulus_size() > ctxt2.coeff_modulus_size()) {
        Plaintext p = encode(1.0, ctxt1);
        mult(p, ctxt1);
    }

    while (ctxt1.coeff_modulus_size() < ctxt2.coeff_modulus_size()) {
        Plaintext p = encode(1.0, ctxt2);
        mult(p, ctxt2);
    }
    //debug
    //if (mode == "debug") {
    //    if (ctxt1.scale() != scales[scales.size() - ctxt1.coeff_modulus_size()]) {
    //        cout << "ctxt1 RESCALING ERROR!!!!" << endl;
    //        cout << "coeff: " << ctxt1.coeff_modulus_size() << " , " << ctxt2.coeff_modulus_size() << endl;
    //        cout << "scale: " << ctxt1.scale() << " , " << ctxt2.scale() << endl;
    //        cout << "---" << endl;
    //    }
    //        
    //    if (ctxt2.scale() != scales[scales.size() - ctxt2.coeff_modulus_size()]) {
    //        cout << "ctxt2 RESCALING ERROR!!!!" << endl;
    //        cout << "coeff: " << ctxt1.coeff_modulus_size() << " , " << ctxt2.coeff_modulus_size() << endl;
    //        cout << "scale: " << ctxt1.scale() << " , " << ctxt2.scale() << endl;
    //        cout << "---" << endl;
    //    }
    //}
}
void ckks_build::scale_equal(Plaintext& ptxt, Ciphertext& ctxt)
{
    Ciphertext x;
    enc->encrypt(ptxt, x);
    scale_equal(x, ctxt);
    dec->decrypt(x, ptxt);
}

/*################ Double Scale Test ################################*/

//multiply plaintext / ciphertext. multiply encoded 1.0 to rescale.
void ckks_build::mul_plain_double(Plaintext& ptxt, Ciphertext& ctxt, Ciphertext& destination)
{
    Plaintext temp = encode(1.0, ctxt);
    eva->multiply_plain(ctxt, ptxt, destination);
    eva->multiply_plain(destination, temp, destination);
    eva->rescale_to_next_inplace(destination);
}
//multiply 3 ciphertexts.
void ckks_build::mul_cipher_double(Ciphertext& ctxt1, Ciphertext& ctxt2, Ciphertext& ctxt3, Ciphertext& destination)
{
    eva->multiply(ctxt1, ctxt2, destination);
    eva->relinearize_inplace(destination, rlk);
    eva->multiply(destination, ctxt3, destination);
    eva->relinearize_inplace(destination, rlk);
    eva->rescale_to_next_inplace(destination);
}

void ckks_build::add_cipher(Ciphertext& ctxt1, Ciphertext& ctxt2, Ciphertext& destination)
{
    eva->add(ctxt1, ctxt2, destination);
}

void ckks_build::temp_d3_doubleScale(vector<double>& poly, Ciphertext& x, Ciphertext& destination)
{
    Plaintext coeff = encode(poly[1], x);
    Plaintext dummy;
    Ciphertext term0, term1;

    mul_plain_double(coeff, x, term0); // ax
    dummy = encode(1.0, term0);
    mul_plain_double(dummy, term0, term0);

    mul_cipher_double(x, x, x, term1); //x^3
    coeff = encode(poly[3], term1);
    mul_plain_double(coeff, term1, term1); // bx^3

    add_cipher(term0, term1, destination); // ax + bx^3
}

void ckks_build::evaluate_function_tripleScale(vector<double>& poly, Ciphertext& x, Ciphertext& destination)
{
    Plaintext coeff;
    Plaintext dummy;
    Ciphertext term0, term1;

    //ax
    coeff = encode(poly[1], x);
    dummy = encode(1.0, x);
    eva->multiply_plain(x, coeff, term0); //1
    eva->multiply_plain(term0, dummy, term0); //2
    eva->multiply_plain(term0, dummy, term0); //3
    eva->rescale_to_next_inplace(term0); // rescale(L-1)

    vector<double> ax = decode_ctxt(term0);
    ax.resize(poly.size());
    printf("ax\n");
    printVector(ax, false, 6);

    //bx^3
    coeff = encode(poly[3], x);
    eva->multiply(x, x, term1); //1
    eva->relinearize_inplace(term1, rlk);
    eva->multiply(x, term1, term1);//2
    eva->relinearize_inplace(term1, rlk);
    eva->multiply_plain(term1, coeff, term1);//3
    eva->rescale_to_next_inplace(term1);// rescale(L-1)

    vector<double> bx = decode_ctxt(term1);
    bx.resize(poly.size());
    printf("bx^3\n");
    printVector(bx, false, 6);

    //ax + bx^3
    eva->add(term0, term1, destination);
}

void ckks_build::evaluate_function_tripleScale_v2(vector<double>& poly, Ciphertext& x, Ciphertext& destination)
{
    Plaintext coeff;
    Plaintext dummy;
    Ciphertext term0, term1;

    int maxScale = 85;
    int curScale;

    //ax
    double a = poly[1] * pow(2.0, 35);
    coeff = encode(a, x);
    eva->multiply_plain(x, coeff, term0);
    term0.scale() = pow(2.0, maxScale);
    eva->rescale_to_next_inplace(term0); // rescale(L-1)

    //bx^3
    double b = poly[3] * pow(2.0, 10);
    coeff = encode(b, x, 1.0);
    term1 = x;
    for (int i = 0; i < 2; i++) {
        eva->multiply(x, term1, term1);
        eva->relinearize_inplace(term1, rlk);
    }
    eva->multiply_plain(term1, coeff, term1);
    term1.scale() = pow(2.0, maxScale);
    eva->rescale_to_next_inplace(term1);// rescale(L-1)

    //ax + bx^3
    eva->add(term0, term1, destination);
}