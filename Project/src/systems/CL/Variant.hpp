#ifndef Variant_HPP
#define Variant_HPP

#include "CL.hpp"


class Variant : public CL
{
    public:
    
    Variant() {};
    Variant(const int WIN, const int MS, const int DS, const int ES) : CL (WIN, MS, DS, ES) {};
    ~Variant() {};
        
    void          keygen    (CLPublickey& clpk,  CLSecretkey& clsk);
    
    // Type 1 Encryption
    CLCiphertext& encrypt      (CLPlaintext& clpt,  CLPublickey& clpk);
    CLPlaintext&  decrypt      (CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk);
    CLPlaintext&  decryptCRT   (CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk);
    CLPlaintext&  hdecrypt     (CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk);
    
    // Type 2 Encryption - Only when fbits <= fmaxbits
    CLCiphertext& encrypt2 (CLPlaintext& clpt,  CLPublickey& clpk);
    CLPlaintext&  decrypt2 (CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk);
    
    // Type 3 Encryption
    CLCiphertext& encrypt3  (CLPlaintext& clpt,  CLPublickey& clpk);
    CLPlaintext&  decrypt3  (CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk);
    CLPlaintext&  hdecrypt3 (CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk);
    
    // Homomorphic functions
    CLCiphertext& evalsum   (CLCiphertext& clct1, CLCiphertext& clct2);
    CLCiphertext& evalscal  (CLCiphertext& clct, const mpz_t alpha);
    CLCiphertext& evalsum3  (CLCiphertext& clct1, CLCiphertext& clct2);
    CLCiphertext& evalscal3 (CLCiphertext& clct, const mpz_t alpha);

}; // end of class

inline void Variant::keygen(CLPublickey& clpk,  CLSecretkey& clsk)
{
    CL::keygen();


	// Set upper bound on size of G
    ZZ_limbs_set(ZZ_UB, B[0]._mp_d, B[0]._mp_size);    


	// Generate sk x in [0, Bf-1] = [0,UB)
    if (ebits == 0)
        ZZ_x = RandomBnd(ZZ_UB);
    else
        ZZ_x = RandomLen_ZZ(ebits);


    // Set G = R^2 of discriminant DK
    G.setD(DK);
    sqr(G, R);
    

    // Compute H = G^x
    qfeG.initialize(G, ZZ_x);
    qfeG.power(H, G, ZZ_x);

    
    // Set exponentiation objects for F and H bases for use in encrypt
    qfeH.initialize(H, e1); // H -- uses qfeH
    
    F.setD(Df);
    qfeF.initialize(F, e1); // F -- uses qfeF

    
    // Set
    clpk.set(N, t);
    clpk.set(UB, f, DK, Df, F, G, H);
    clsk.set(ZZ_x);
}


// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------
// ----------------------------------- TYPE 1 -----------------------------------
// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------

inline CLCiphertext& Variant::encrypt(CLPlaintext& clpt, CLPublickey& clpk)
{
    // Get pk, m
    clpk.get(UB, f, Df, G, H);
    clpk.get(N, t);
    clpt.get(m, ZZ_m);

    // Select r from [0, Bcon-1]
    if (ebits == 0)
        e1 = RandomBnd(ZZ_UB);
    else
        e1 = RandomLen_ZZ(ebits);

	
    // Compute C1 = G^r ; has D = DK
    G.setD(DK);
    qfeG.power(T1, G, e1);
    
    // Compute C2 = F^m * psi(H^r)
    qfeH.power(T4, H, e1);
    liftQF(T5, T4, DK, Df, f); // phi_inv(H^r)
    
    
    if (t == 1)
    {
        qfe.initialize(T5, ZZ_f);
        qfe.power(T3, T5, ZZ_f); // psi(H^r)
        
        qfe.CLpower(T2, f, Df, m); // F^m
        mul(T3, T3, T2);
    }
    else
    {
        qfed.initialize(F, T5, ZZ_m, ZZ_f);
        qfed.power(T3, F, T5, ZZ_m, ZZ_f);
    }
    
    
    // Set ciphertext
    ctt.set(T1, T3);
    
    // Return ciphertext
    return ctt;
}


inline CLPlaintext& Variant::decrypt(CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk)
{
    // Get pk, sk, ciphertexts
    clpk.get(f);
    clsk.get(ZZ_x);
    clct.get(C1, C2);
    
    // Compute M = C2 * psi(C1^-x)
    C1.setD(DK);
    qfe.initialize(C1, ZZ_x);
    qfe.power(M, C1, ZZ_x);
    inv(M, M); // C1^-x
    liftQF(T1, M, DK, Df, f);
    qfe.initialize(T1, ZZ_f);
    qfe.power(M, T1, ZZ_f); // psi(C1^-x)
    mul(M, M, C2); // C2 * psi(C1^-x)
    
    
    // Get decrypted m from M
    CL::solve(dm, M, f);
    
    // Set decrypted result
    ptt.set(dm);
    
    // Return decrypted plaintext
    return ptt;
}


inline CLPlaintext& Variant::decryptCRT(CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk)
{
    // Get pk, sk, ciphertexts
    clpk.get(f);
    clsk.get(ZZ_x);
    clct.get(C1, C2);
    
    // Compute M = C2 * psi(C1^-x)
    C1.setD(DK);
    qfe.initialize(C1, ZZ_x);
    qfe.power(M, C1, ZZ_x);
    inv(M, M); // C1^-x
    liftQF(T1, M, DK, Df, f);
    qfe.initialize(T1, ZZ_f);
    qfe.power(M, T1, ZZ_f); // psi(C1^-x)
    mul(M, M, C2);
    
    // Get L(m) part of M
    M.getb(temp1);
    mpz_divexact(temp1, temp1, f);
    
    // Compute L(m)^-1 mod pi^t
    for (int i = 0; i < N; i++)
        mpz_invert(var[i], temp1, pit[i]);
    
    // Apply CRT to get m
    CL::CRTl(dm, N, var, pit);
    
    // Set decrypted result
    ptt.set(dm);
    
    // Return decrypted plaintext
    return ptt;
}


inline CLPlaintext& Variant::hdecrypt(CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk)
{
    // Get pk, sk, ciphertexts
    clpk.get(f);
    clsk.get(ZZ_x);
    clct.get(C1, C2);
    
    // Compute M = C2 * psi(C1^-x)
    C1.setD(DK);
    qfe.initialize(C1, ZZ_x);
    qfe.power(M, C1, ZZ_x);
    inv(M, M); // C1^-x
    liftQF(T1, M, DK, Df, f);
    qfe.initialize(T1, ZZ_f);
    qfe.power(M, T1, ZZ_f); // psi(C1^-x)
    mul(M, M, C2); // C2 * psi(C1^-x)
    
    // Get decrypted m from M
    switch(N)
    {
        case (1):
            CL::solve(dm, F, M, f, p, t);  // Prime power PH
            break;
            
        default:
            CL::solve(dm, F, M, f, pit, pN, N, t);  // PH
            break;
    }
    
    // Set decrypted result
    ptt.set(dm);
    
    // Return decrypted plaintext
    return ptt;
}


// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------
// ----------------------------------- TYPE 2 -----------------------------------
// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------


inline CLCiphertext& Variant::encrypt2(CLPlaintext& clpt, CLPublickey& clpk)
{
    // Get pk, m
    clpk.get(UB, f, Df, G, H, pN, N);
    clpk.get(N,t);
    clpt.get(m);
    
    // Compute m1 = m mod p1^t and m2 = m mod p2^t
    for (int i = 0; i < N; i++)
        mpz_mod(mN[i], m, pit[i]);
    
    
    // Select r from [0, Bcon-1]
    if (ebits == 0)
        e1 = RandomBnd(ZZ_UB);
    else
        e1 = RandomLen_ZZ(ebits);
    
    
    // Compute C1 = G^r
    G.setD(DK);
    qfeG.power(T1, G, e1);
    
    
    // Compute C2 = F1^m1 * ... * FN^mN * psi(H^r)
    qfeH.power(T4, H, e1);
    liftQF(T5, T4, DK, Df, f);
    qfe.initialize(T5, ZZ_f);
    qfe.power(T2, T5, ZZ_f); // psi(H^r)
    
    if (fbits <= fmaxbits)
    {
        for (int i = 0; i < N; i++)
        {
            qfe.CLpower(TN[i], pit[i], Df, mN[i]);
            mul(T2, T2, TN[i]);
        }
    }
    
    else
    {
        for (int i = 0; i < N; i++)
        {
            ZZ_limbs_set(e1, mN[i][0]._mp_d, mN[i][0]._mp_size);
            qfe.initialize(FN[i], e1);
            qfe.power(TN[i], FN[i], e1);
            mul(T2, T2, TN[i]);
        }
    }
    
    // Set ciphertext
    ctt.set(T1, T2);
    
    // Return ciphertext
    return ctt;
}


inline CLPlaintext& Variant::decrypt2(CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk)
{
    // Get pk, sk, ciphertexts
    clpk.get(f, DK, Df, pN, N);
    clsk.get(ZZ_x);
    clct.get(C1, C2);
    
    // Compute M = C2/psi(C1^x) = (N^2, X)
    C1.setD(DK);
    qfe.initialize(C1, ZZ_x);
    qfe.power(M, C1, ZZ_x);
    inv(M, M); // C1^-x
    liftQF(T1, M, DK, Df, f);
    qfe.initialize(T1, ZZ_f);
    qfe.power(M, T1, ZZ_f); // psi(C1^-x)
    mul(M, M, C2);
    
    // Get b part of M = (N^2, X)
    M.getb(temp1);
    
    // Compute (X/pi^t)^-1 mod pi^t
    for (int i = 0; i < N; i++)
    {
        mpz_divexact(temp2, temp1, pit[i]);
        mpz_invert(var[i], temp2, pit[i]);
    }
    
    // Apply CRT to get m
    CL::CRTl(dm, N, var, pit);
    ptt.set(dm);
    
    // Return decrypted plaintext
    return ptt;
}


// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------
// ----------------------------------- TYPE 3 -----------------------------------
// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------


inline CLCiphertext& Variant::encrypt3(CLPlaintext& clpt, CLPublickey& clpk)
{
    // Get pk, m
    clpk.get(UB, f, DK, Df, F, G, H);
    clpk.get(N,t);
    clpt.get(m);
    
    // Compute mN = m mod pN^tn
    for (int i = 0; i < N; i++)
        mpz_mod(mN[i], m, pit[i]);
    
    
    // Select r from [0, Bcon-1]
    if (ebits == 0)
        e1 = RandomBnd(ZZ_UB);
    else
        e1 = RandomLen_ZZ(ebits);
    
    
    // Compute C1 = G^r
    G.setD(DK);
    qfeG.power(T1, G, e1);
    
    
    // Compute CN = FN^mN * psi(H^r)
    qfeH.power(T4, H, e1);
    liftQF(T5, T4, DK, Df, f);
    qfe.initialize(T5, ZZ_f);
    qfe.power(T2, T5, ZZ_f); // psi(H^r)
    
    
    if (t == 1)
    {
        for (int i = 0; i < N; i++)
        {
            qfe.CLpower(CN[i], pit[i], Df, mN[i]);
            mul(CN[i], CN[i], T2);
        }
    }
    
    else
    {
        for (int i = 0; i < N; i++)
        {
            ZZ_limbs_set(e2, mN[i]->_mp_d, mN[i]->_mp_size);
            qfeFN[i].power(CN[i], FN[i], e2);
            mul(CN[i], CN[i], T2);
        }
    }
    
    // Set ciphertext
    ctt.set(T1, CN, N);
    
    
    // Return ciphertext
    return ctt;
}



inline CLPlaintext& Variant::decrypt3(CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk)
{
    // Get pk, sk, ciphertexts
    clpk.get(F, f, pN, N, t);
    clsk.get(ZZ_x);
    clct.get(C1, CN, N);
    
    // Compute psi(C1^-x)
    C1.setD(DK);
    qfe.initialize(C1, ZZ_x);
    qfe.power(CIX, C1, ZZ_x);
    inv(CIX, CIX);
    liftQF(T1, CIX, DK, Df, f);
    qfe.initialize(T1, ZZ_f);
    qfe.power(CIX, T1, ZZ_f); // CIX = psi(C1^-x)
    
    // Compute mN from Mn = CN/psi(C1^x)
    for (int i = 0; i < N; i++)
    {
        mul(M, CN[i], CIX);
        CL::solve(var[i], M, pit[i]);
    }
    
    // CRT
    CL::CRTl(dm, N, var, pit);
    ptt.set(dm);
    
    
    // Return decrypted plaintext
    return ptt;
}

inline CLPlaintext& Variant::hdecrypt3(CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk)
{
    // Get pk, sk, ciphertexts
    clpk.get(F, f, pN, N, t);
    clsk.get(ZZ_x);
    clct.get(C1, CN, N);
    
    // Compute psi(C1^-x)
    C1.setD(DK);
    qfe.initialize(C1, ZZ_x);
    qfe.power(CIX, C1, ZZ_x);
    inv(CIX, CIX);
    liftQF(T1, CIX, DK, Df, f);
    qfe.initialize(T1, ZZ_f);
    qfe.power(CIX, T1, ZZ_f); // CIX = psi(C1^-x)
    
    // Compute mN from Mn = CN/psi(C1^x)
    for (int i = 0; i < N; i++)
    {
        mul(M, CN[i], CIX);
        CL::solve(var[i], FN[i], M, pit[i], pN[i], t);
    }
    
    // CRT
    CL::CRTl(dm, N, var, pit);
    ptt.set(dm);
    
    
    // Return decrypted plaintext
    return ptt;
}


// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------
// ----------------- HOMOMORPHIC PROPERTY VALIDATION FUNCTIONS ------------------
// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------


inline CLCiphertext& Variant::evalsum(CLCiphertext& clct1, CLCiphertext& clct2)
{
    clct1.get(T1, T2);
    clct2.get(T3, T4);

    T6.setD(Df);
    mul(T6, T2, T4);
    
    T5.setD(DK);
    mul(T5, T1, T3);
    
    // Select r from [0, Bcon-1]
    if (ebits == 0)
        e1 = RandomBnd(ZZ_UB);
    else
        e1 = RandomLen_ZZ(ebits);
        
    // Compute G^r
    qfeG.power(T1, G, e1);

    mul(T5, T5, T1);
    
    // Compute psi(H^r)
    qfeH.power(T2, H, e1);
    liftQF(T4, T2, DK, Df, f);
    qfe.initialize(T4, ZZ_f);
    qfe.power(T2, T4, ZZ_f);
        
    mul(T6, T6, T2);
    
    // Set ciphertext
    ctt.set(T5, T6);
        
    // Return ciphertext
    return ctt;
}


inline CLCiphertext& Variant::evalscal(CLCiphertext& clct, const mpz_t alpha)
{
    clct.get(T1, T2);
    ZZ_limbs_set(e1, alpha[0]._mp_d, alpha[0]._mp_size);
    
    T1.setD(DK);
    qfe.initialize(T1, e1);
    qfe.power(T3, T1, e1);
    
    T2.setD(Df);
    qfe2.initialize(T2, e1);
    qfe2.power(T4, T2, e1);
    
    // Set ciphertext
    ctt.set(T3, T4);
    
    // Return ciphertext
    return ctt;
}



inline CLCiphertext& Variant::evalsum3(CLCiphertext& clct1, CLCiphertext& clct2)
{
    clct1.get(T1, TN1, N);
    clct2.get(T3, TN2, N);
    
    T5.setD(DK);
    mul(T5, T1, T3);

    // Select r from [0, Bcon-1]
    if (ebits == 0)
        e1 = RandomBnd(ZZ_UB);
    else
        e1 = RandomLen_ZZ(ebits);

    // Compute G^r
    qfeG.power(T1, G, e1);
    mul(T5, T5, T1);

    qfeH.power(T2, H, e1);
    
    // Compute psi(H^r)
    T4.setD(Df);
    liftQF(T4, T2, DK, Df, f);
    qfe.initialize(T4, ZZ_f);
    qfe.power(T6, T4, ZZ_f);
    for (int i = 0; i < N; i++)
    {
        mul(TN3[i], TN1[i], TN2[i]);
    }
    
    for (int i = 0; i < N; i++)
    {
        mul(TN3[i], TN3[i], T6);
    }
    
    // Set ciphertext
    ctt.set(T5, TN3, N);
    
    // Return ciphertext
    return ctt;
}


inline CLCiphertext& Variant::evalscal3(CLCiphertext& clct, const mpz_t alpha)
{
    clct.get(T1, TN1, N);
    ZZ_limbs_set(e1, alpha[0]._mp_d, alpha[0]._mp_size);
    
    T3.setD(DK);
    qfe.initialize(T1, e1);
    qfe.power(T3, T1, e1);
    T3.setD(Df);
    for (int i = 0; i < N; i++)
    {
        qfeTN[i].initialize(TN1[i], e1);
        qfeTN[i].power(TN2[i], TN1[i], e1); // FN^mN
    }
    
    ctt.set(T3, TN2, N);
    return ctt;
}

#endif // Variant_HPP

