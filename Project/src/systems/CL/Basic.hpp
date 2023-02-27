#ifndef Basic_HPP
#define Basic_HPP

#include "CL.hpp"


class Basic : public CL
{
    public:
    
    Basic() {};
    Basic (const int WIN, const int MS, const int DS, const int ES) : CL (WIN, MS, DS, ES) {};
    ~Basic() {};

    
    void          keygen      (CLPublickey& clpk,  CLSecretkey& clsk);
    
    //
    CLCiphertext& encrypt      (CLPlaintext& clpt,  CLPublickey& clpk);
    CLPlaintext&  decrypt      (CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk);
    CLPlaintext&  decryptCRT   (CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk);
    CLPlaintext&  hdecrypt     (CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk);
    
    //
    CLCiphertext& encrypt2 (CLPlaintext& clpt,  CLPublickey& clpk);
    CLPlaintext&  decrypt2 (CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk);
    
    //
    CLCiphertext& encrypt3  (CLPlaintext& clpt,  CLPublickey& clpk);
    CLPlaintext&  decrypt3  (CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk);
    CLPlaintext&  hdecrypt3 (CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk);
    
    // Homomorphic functions
    CLCiphertext& evalsum   (CLCiphertext& clct1, CLCiphertext& clct2);
    CLCiphertext& evalscal  (CLCiphertext& clct, const mpz_t alpha);
    CLCiphertext& evalsum3  (CLCiphertext& clct1, CLCiphertext& clct2);
    CLCiphertext& evalscal3 (CLCiphertext& clct, const mpz_t alpha);

}; // end of class

inline void Basic::keygen(CLPublickey& clpk, CLSecretkey& clsk)
{
    CL::keygen();


	// Set upper bound on size of G
    mpz_mul(UB, B, f);
    ZZ_limbs_set(ZZ_UB, UB[0]._mp_d, UB[0]._mp_size);

    
	// Generate sk x in [0, Bf-1] = [0,UB)
    if (ebits == 0)
        ZZ_x = RandomBnd(ZZ_UB);
    else
        ZZ_x = RandomLen_ZZ(ebits);


    // Set G = phi_inv(R^2)^f * (F^k)
    G.setD(Df);
    ZZ_limbs_set(e1, k[0]._mp_d, k[0]._mp_size);
    qfed.initialize(psi_RR, F, ZZ_f, e1);
    qfed.power(G, psi_RR, F, ZZ_f, e1);

    
    // Compute H = G^x
    qfeG.initialize(G, ZZ_x);
    qfeG.power(H, G, ZZ_x);


    // Set exponentiation objects for F and H bases for use in encrypt
    qfed.initialize(F, H, e1, e1); // Double exp.
    qfeF.initialize(F, e1); // F -- uses qfeF
    qfeH.initialize(H, e1); // H -- uses qfeH

    
    // Set public and secret key
    clpk.set(N, t);
    clpk.set(UB, f, DK, Df, F, G, H);
    clsk.set(ZZ_x);
}



inline CLCiphertext& Basic::encrypt(CLPlaintext& clpt, CLPublickey& clpk)
{
    // Get pk, m
    clpk.get(UB, f, Df, G, H);
    clpk.get(N, t);
    clpt.get(m, ZZ_m);
    
    // Select r from [0, Bf-1]
    if (ebits == 0)
        e1 = RandomBnd(ZZ_UB);
    else
        e1 = RandomLen_ZZ(ebits);
    
    // Compute C1 = G^r
    qfeG.power(T1, G, e1);
    
    // Compute C2 = F^m * H^r
    if (t == 1)
    {
        qfeH.power(T3, H, e1); // H^r
        
        qfe.CLpower(T2, f, Df, m); // F^m
        mul(T3, T3, T2); // F^m * H^r
    }
    else
        qfed.power(T3, F, H, ZZ_m, e1);
    
    
    // Set ciphertext C1, C2
    ctt.set(T1, T3);
    
    
    // Return ciphertext
    return ctt;
}


inline CLPlaintext& Basic::decrypt(CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk)
{
    // Get pk, sk, ciphertexts
    clpk.get(f);
    clsk.get(ZZ_x);
    clct.get(C1, C2);
    
    // Compute M = C2/C1^x
    qfe.initialize(C1, ZZ_x);
    qfe.power(CIX, C1, ZZ_x);
    inv(CIX, CIX);
    mul(M, C2, CIX);
    
    // Get decrypted m from M - Modular inverse decryption
    CL::solve(dm, M, f);
    
    // Set decrypted result
    ptt.set(dm);
    
    // Return decrypted plaintext
    return ptt;
}


inline CLPlaintext& Basic::decryptCRT(CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk)
{
    // Get pk, sk, ciphertexts
    clpk.get(f);
    clsk.get(ZZ_x);
    clct.get(C1, C2);
    
    // Compute M = C2/C1^x
    qfe.initialize(C1, ZZ_x);
    qfe.power(CIX, C1, ZZ_x);
    inv(CIX, CIX);
    mul(M, C2, CIX);
    
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


inline CLPlaintext& Basic::hdecrypt(CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk)
{
    // Get pk, sk, ciphertexts
    clpk.get(f);
    clsk.get(ZZ_x);
    clct.get(C1, C2);
    
    // Compute M = C2/C1^x
    qfe.initialize(C1, ZZ_x);
    qfe.power(CIX, C1, ZZ_x);
    inv(CIX, CIX);
    mul(M, C2, CIX);
    
    switch(N)
    {
        case (1): // Prime power PH
            CL::solve(dm, F, M, f, p, t);
            break;
            
        default: // Generic PH
            CL::solve(dm, F, M, f, pit, pN, N, t);
            break;
    }
    
    // Set decrypted result
    ptt.set(dm);
    
    
    // Return decrypted plaintext
    return ptt;
}




inline CLCiphertext& Basic::encrypt2(CLPlaintext& clpt, CLPublickey& clpk)
{
    // Get pk, m
    clpk.get(UB, f, Df, G, H, pN, N);
    clpk.get(N, t);
    clpt.get(m);
    
    // Compute m1 = m mod p1^t and m2 = m mod p2^t
    for (int i = 0; i < N; i++)
        mpz_mod(mN[i], m, pit[i]);
    
    // Select r from [0, Bf-1]
    if (ebits == 0)
        e1 = RandomBnd(ZZ_UB);
    else
        e1 = RandomLen_ZZ(ebits);
    
    // Compute C1 = G^r
    qfeG.power(T1, G, e1);
    
    // Compute C2 = F1^m1 * ... * FN^mN * H^r
    qfeH.power(T4, H, e1);
    if (fbits <= fmaxbits)
    {
        for (int i = 0; i < N; i++)
        {
            qfe.CLpower(TN[i], pit[i], Df, mN[i]);
            mul(T4, T4, TN[i]);
        }
    }
    else
    {
        for (int i = 0; i < N; i++)
        {
            ZZ_limbs_set(e1, mN[i][0]._mp_d, mN[i][0]._mp_size);
            qfe.initialize(FN[i], e1);
            qfe.power(TN[i], FN[i], e1);
            mul(T4, T4, TN[i]);
        }
    }
    
    // Set ciphertext
    ctt.set(T1, T4);
    
    // Return ciphertext
    return ctt;
}


inline CLPlaintext& Basic::decrypt2(CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk)
{
    // Get pk, sk, ciphertexts
    clpk.get(f, pN, N);
    clsk.get(ZZ_x);
    clct.get(C1, C2);
    
    // Compute M = C2/C1^x = (N^2, X)
    qfe.initialize(C1, ZZ_x);
    qfe.power(CIX, C1, ZZ_x);
    inv(CIX, CIX);
    mul(M, C2, CIX);
    
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


inline CLCiphertext& Basic::encrypt3(CLPlaintext& clpt, CLPublickey& clpk)
{
    // Get pk, m
    clpk.get(UB, f, F, G, H);
    clpk.get(N,t);
    clpt.get(m);
    
    // Compute mN = m mod pN^tn
    for (int i = 0; i < N; i++)
        mpz_mod(mN[i], m, pit[i]);
    
    
    // Select r from [0, Bf-1]
    if (ebits == 0)
        e1 = RandomBnd(ZZ_UB);
    else
        e1 = RandomLen_ZZ(ebits);
    
    
    // Compute C1 = G^r
    qfeG.power(T1, G, e1);
    
    
    // Compute Cn = FN^mN * H^r
    qfeH.power(T2, H, e1);
    
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
            ZZ_limbs_set(e2, mN[i][0]._mp_d, mN[i][0]._mp_size);
            qfeFN[i].power(CN[i], FN[i], e2);
            mul(CN[i], CN[i], T2);
        }
    }
    
    // Set ciphertext
    ctt.set(T1, CN, N);
    
    
    // Return ciphertext
    return ctt;
}


inline CLPlaintext& Basic::decrypt3(CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk)
{
    // Get pk, sk, ciphertexts
    clpk.get(F, f, pN, N, t);
    clsk.get(ZZ_x);
    clct.get(C1, CN, N);
    
    // Compute C1^-x
    qfe.initialize(C1, ZZ_x);
    qfe.power(CIX, C1, ZZ_x);
    inv(CIX, CIX);
    
    // Compute mN from MN = Cn/C1^x
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


inline CLPlaintext& Basic::hdecrypt3(CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk)
{
    // Get pk, sk, ciphertexts
    clpk.get(F, f, pN, N, t);
    clsk.get(ZZ_x);
    clct.get(C1, CN, N);
    
    // Compute C1^-x
    qfe.initialize(C1, ZZ_x);
    qfe.power(CIX, C1, ZZ_x);
    inv(CIX, CIX);
    
    // Compute mN from MN = Cn/C1^x
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


inline CLCiphertext& Basic::evalsum(CLCiphertext& clct1, CLCiphertext& clct2)
{
    clct1.get(T1, T2);
    clct2.get(T3, T4);
    
    //
    mul(T5, T1, T3);
    mul(T6, T2, T4);
    
    // Select r from [0, Bf-1]
    if (ebits == 0)
        e1 = RandomBnd(ZZ_UB);
    else
        e1 = RandomLen_ZZ(ebits);
        
    // Compute G^r
    qfeG.power(T1, G, e1);
    
    // Compute H^r
    qfeH.power(T2, H, e1);
    
    //
    mul(T5, T5, T1);
    mul(T6, T6, T2);
    
    // Set ciphertext
    ctt.set(T5, T6);
        
    // Return ciphertext
    return ctt;
}


inline CLCiphertext& Basic::evalscal(CLCiphertext& clct, const mpz_t alpha)
{
    clct.get(T1, T2);
    ZZ_limbs_set(e1, alpha[0]._mp_d, alpha[0]._mp_size);
    
    qfe .initialize(T1, e1);
    qfe2.initialize(T2, e1);
    
    qfe .power(T3, T1, e1);
    qfe2.power(T4, T2, e1);
    
    // Set ciphertext
    ctt.set(T3, T4);
    
    // Return ciphertext
    return ctt;
}


inline CLCiphertext& Basic::evalsum3(CLCiphertext& clct1, CLCiphertext& clct2)
{
    
    clct1.get(T1, TN1, N);
    clct2.get(T3, TN2, N);
    
    mul(T5, T1, T3);
    for (int i = 0; i < N; i++)
    {
        mul(TN3[i], TN1[i], TN2[i]);
    }
    
    // Select r from [0, Bf-1]
    if (ebits == 0)
        e1 = RandomBnd(ZZ_UB);
    else
        e1 = RandomLen_ZZ(ebits);
    
    // Compute G^r
    qfeG.power(T1, G, e1);
    
    // Compute H^r
    qfeH.power(T2, H, e1);
    
    //
    mul(T5, T5, T1);
    
    for (int i = 0; i < N; i++)
    {
        mul(TN3[i], TN3[i], T2);
    }
    
    // Set ciphertext
    ctt.set(T5, TN3, N);
    
    // Return ciphertext
    return ctt;
}


inline CLCiphertext& Basic::evalscal3(CLCiphertext& clct, const mpz_t alpha)
{
    clct.get(T1, TN1, N);
    ZZ_limbs_set(e1, alpha[0]._mp_d, alpha[0]._mp_size);
    
    qfe.initialize(T1, e1);
    qfe.power(T3, T1, e1);
    for (int i = 0; i < N; i++)
    {
        qfeTN[i].initialize(TN1[i], e1);
        qfeTN[i].power(TN2[i], TN1[i], e1); // FN^mN
    }
    
    ctt.set(T3, TN2, N);
    return ctt;
}

#endif // Basic_HPP
