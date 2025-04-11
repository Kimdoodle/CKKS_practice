#include "header/SEAL_MM.h"

//Generator
ckks_build::ckks_build(double scale, size_t pmd)
{
    this->scale = scale;
    this->poly_modulus_degree = pmd;

    parms = make_unique<EncryptionParameters>(scheme_type::ckks);
    parms->set_poly_modulus_degree(pmd);
    
    // set modulus chain.
    modulus = { 60, 60, 60 , 60 , 60, 60 };

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
}

void ckks_build::createGaloisKeys(int start, int end)
{
    vector<int> list;
    for (int i = start; i <= end; i++)
        list.push_back(i);
    keygen->create_galois_keys(list, glk);
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

Plaintext ckks_build::encode(const vector<double>& input, Ciphertext& ctxt)
{
    Plaintext plain;
    encoder->encode(input, ctxt.parms_id(), ctxt.scale(), plain);
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

//equal two ciphertext's modulus level.
//void ckks_build::modulus_equal(Ciphertext& ctxt1, Ciphertext& ctxt2)
//{
//    auto level1 = ctxt1.coeff_modulus_size();
//    auto level2 = ctxt2.coeff_modulus_size();
//
//    if (level1 < level2) {
//        modulus_switch(ctxt2, ctxt1.parms_id());
//    }
//    else if (level1 > level2) {
//        modulus_switch(ctxt1, ctxt2.parms_id());
//    }
//}

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
//void ckks_build::mul_plain_double(Plaintext& ptxt, Ciphertext& ctxt, Ciphertext& destination)
//{
//    Plaintext temp = encode(1.0, ctxt);
//    eva->multiply_plain(ctxt, ptxt, destination);
//    eva->multiply_plain(destination, temp, destination);
//    eva->rescale_to_next_inplace(destination);
//}
////multiply 3 ciphertexts.
//void ckks_build::mul_cipher_double(Ciphertext& ctxt1, Ciphertext& ctxt2, Ciphertext& ctxt3, Ciphertext& destination)
//{
//    eva->multiply(ctxt1, ctxt2, destination);
//    eva->relinearize_inplace(destination, rlk);
//    eva->multiply(destination, ctxt3, destination);
//    eva->relinearize_inplace(destination, rlk);
//    eva->rescale_to_next_inplace(destination);
//}
//
//void ckks_build::add_cipher(Ciphertext& ctxt1, Ciphertext& ctxt2, Ciphertext& destination)
//{
//    eva->add(ctxt1, ctxt2, destination);
//}
//
//void ckks_build::temp_d3_doubleScale(vector<double>& poly, Ciphertext& x, Ciphertext& destination)
//{
//    Plaintext coeff = encode(poly[1], x);
//    Plaintext dummy;
//    Ciphertext term0, term1;
//
//    mul_plain_double(coeff, x, term0); // ax
//    dummy = encode(1.0, term0);
//    mul_plain_double(dummy, term0, term0);
//
//    mul_cipher_double(x, x, x, term1); //x^3
//    coeff = encode(poly[3], term1);
//    mul_plain_double(coeff, term1, term1); // bx^3
//
//    add_cipher(term0, term1, destination); // ax + bx^3
//}
//
//void ckks_build::evaluate_function_tripleScale(vector<double>& poly, Ciphertext& x, Ciphertext& destination)
//{
//    Plaintext coeff;
//    Plaintext dummy;
//    Ciphertext term0, term1;
//
//    //ax
//    coeff = encode(poly[1], x);
//    dummy = encode(1.0, x);
//    eva->multiply_plain(x, coeff, term0); //1
//    eva->multiply_plain(term0, dummy, term0); //2
//    eva->multiply_plain(term0, dummy, term0); //3
//    eva->rescale_to_next_inplace(term0); // rescale(L-1)
//
//    vector<double> ax = decode_ctxt(term0);
//    ax.resize(poly.size());
//    printf("ax\n");
//    printVector(ax, false, 6);
//
//    //bx^3
//    coeff = encode(poly[3], x);
//    eva->multiply(x, x, term1); //1
//    eva->relinearize_inplace(term1, rlk);
//    eva->multiply(x, term1, term1);//2
//    eva->relinearize_inplace(term1, rlk);
//    eva->multiply_plain(term1, coeff, term1);//3
//    eva->rescale_to_next_inplace(term1);// rescale(L-1)
//
//    vector<double> bx = decode_ctxt(term1);
//    bx.resize(poly.size());
//    printf("bx^3\n");
//    printVector(bx, false, 6);
//
//    //ax + bx^3
//    eva->add(term0, term1, destination);
//}

void ckks_build::evaluate_function_tripleScale_v2(vector<double>& poly, Ciphertext& x, Ciphertext& destination)
{
    Plaintext coeff;
    Plaintext dummy;
    Ciphertext term0, term1;

    int maxScale = 85;
    int curScale;

    /*
        1 rescaling after 3 multiplication.
        Scale size = 2^25
        modulus size = 60 bits
        ex) ax
            coeff a size = 1.0(scale) * 2^60
            ciphertext x size = 2^25(scale)
            a * x  size = 2^85

        ex2) bx^3
            coeff b size = 1.0(scale) * 2^10
            ciphertext x^3 size = (2^25)^3 = 2^75(scale)
            bx^3 size = 2^85(scale)
    */
    //ax
    double a = poly[1] * pow(2.0, 60); 
    coeff = encode(a, x, 1.0);
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

/* ######################################################## */

vector<Plaintext> ckks_build::encode_matrix(vector<vector<double>> U, Ciphertext& x)
{
    vector<Plaintext> temp;
    for (vector<double> row : U)
    {
        temp.push_back(encode(row, x));
    }
    return temp;
}

Ciphertext ckks_build::rotate(Ciphertext &x, int k)
{
    Ciphertext res;
    eva->rotate_vector(x, k, glk, res);
    return res;
}

Ciphertext ckks_build::linear_transform_sub(Ciphertext& x, Plaintext& coeffs, int rot, int xLength)
{
    Ciphertext rotated_x, temp1, temp2, result;
    rotated_x = rotate(x, rot);
    eva->multiply_plain(rotated_x, coeffs, temp1);
    rotated_x = rotate(x, rot - xLength);
    eva->multiply_plain(rotated_x, coeffs, temp2);
    eva->add(temp1, temp2, result);

    //cout << "----------------------- " << rot << " ---------------------- - " << endl;
    //debug_vector(decode_ctxt(x), 9, "x<-" + to_string(rot));
    //debug_vector(decode_ctxt(result), 9, "ctprime");

    return result;
}

Ciphertext ckks_build::linear_transform(Ciphertext& x, vector<vector<double>>& dVecs, int xLength)
{
    Ciphertext ctprime, temp, rotated_x;
    rotated_x = x;
    //1이 포함된 대각벡터만 추출
    vector<vector<double>> valid_dVecs;
    vector<int> valid_rotations;

    for (int i = 0; i < dVecs.size(); ++i) {
        const auto& row = dVecs[i];
        bool isValid = std::any_of(row.begin(), row.end(), [](double val) {
            return val != 0.0;
            });
        if (isValid) {
            valid_dVecs.push_back(row);
            valid_rotations.push_back(i);
        }
    }

    vector<Plaintext> encoded_dVecs = encode_matrix(valid_dVecs, x);
    ctprime = linear_transform_sub(x, encoded_dVecs[0], valid_rotations[0], xLength);

    for (int i = 1; i < encoded_dVecs.size(); i++)
    {
        temp = linear_transform_sub(x, encoded_dVecs[i], valid_rotations[i], xLength);
        eva->add(temp, ctprime, ctprime);
    }
    return ctprime;
}

Ciphertext ckks_build::matrix_multiplication(Ciphertext& A, Ciphertext& B, int originalD)
{
    createGaloisKeys(-originalD, originalD);

    vector<vector<double>> matrix, d_matrix, pad_d_matrix;
    vector<Plaintext> encoded_matrix;
    int sizeU = int(sqrt(poly_modulus_degree / 2));
    int d = int(sqrt(originalD));

    // Step 1-1
    matrix = make_matrix(d, "sigma");
    Ciphertext ctA0 = linear_transform(A, matrix, originalD); // scale+1
    eva->rescale_to_next_inplace(ctA0);
    //debug_printMatrix(ctA0, originalD, "ctA0");

    //Step 1-2
    matrix = make_matrix(d, "tau");
    Ciphertext ctB0 = linear_transform(B, matrix, originalD); // scale+1
    eva->rescale_to_next_inplace(ctB0);
    //debug_printMatrix(ctB0, originalD, "ctB0");

    //Step 2
    vector<Ciphertext> ctAk_vector, ctBk_vector;
    Ciphertext ctAk, ctBk;
    for (int k = 0; k < originalD; k++)
    {
        matrix = make_matrix(d, "phi" + to_string(k));
        ctAk = linear_transform(ctA0, matrix, originalD);
        eva->rescale_to_next_inplace(ctAk);
        ctAk_vector.push_back(ctAk);

        matrix = make_matrix(d, "psi" + to_string(k));
        ctBk = linear_transform(ctB0, matrix, originalD);
        eva->rescale_to_next_inplace(ctBk);
        ctBk_vector.push_back(ctBk);

        //debug_printMatrix(ctAk, originalD, "ctA" + to_string(k));
        //debug_printMatrix(ctBk, originalD, "ctB" + to_string(k));
    }
    //eva->rescale_to_next_inplace(ctA0);
    //eva->rescale_to_next_inplace(ctB0);
    //ctAk_vector.insert(ctAk_vector.begin(), ctA0);
    //ctBk_vector.insert(ctBk_vector.begin(), ctB0);

    //Step3
   Ciphertext ctAB, temp;
    eva->multiply(ctAk_vector[0], ctBk_vector[0], ctAB);
    eva->relinearize_inplace(ctAB, rlk);

    for (int k = 1; k < originalD; k++)
    {
        eva->multiply(ctAk_vector[k], ctBk_vector[k], temp);
        eva->relinearize_inplace(temp, rlk);
        eva->add(ctAB, temp, ctAB);

        //debug_vector(decode_ctxt(temp), 9, "ctA" + to_string(k));
        //debug_vector(decode_ctxt(ctAB), 9, "ctAB");
    }
    eva->rescale_to_next_inplace(ctAB);
    return ctAB;
}

void ckks_build::debug_printMatrix(Ciphertext& A, int originalD, string title)
{
    vector<double> decA = decode_ctxt(A);
    vector<vector<double>> unflat_decA = unflatten_matrix(decA, int(sqrt(decA.size())));
    vector<vector<double>> final_decA = ipad_matrix(unflat_decA, originalD);
    cout << title << endl;
    printMatrix(final_decA);
}

void ckks_build::debug_vector(vector<double> A, int size, string title)
{
    cout << title << endl;
    vector<double> dectest = A;
    dectest.resize(size);
    printVector(dectest, false, 0);
}