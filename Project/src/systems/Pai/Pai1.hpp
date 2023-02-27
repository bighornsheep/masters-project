/**
 *
 * AUTHOR - Parthasarathi Das
 * BRIEF  - A class for the Paillier (variant) cryptosystem
 */

#ifndef Pai1_HPP
#define Pai1_HPP

#include "Pai.hpp"

class Pai1 : public Pai
{
    public:
     Pai1(): Pai() {};
    ~Pai1() {};
    
    void           keygen     (PaiPublickey& paipk, PaiSecretkey& paisk);
    PaiCiphertext& encrypt    (PaiPlaintext& paipt, PaiPublickey& paipk);
    PaiPlaintext&  decrypt    (PaiCiphertext& paict, PaiPublickey& paipk, PaiSecretkey& paisk);
    PaiPlaintext&  decryptCRT (PaiCiphertext& paict, PaiPublickey& paipk, PaiSecretkey& paisk);
    
}; // end of class


// Out of line class member definitions
inline void Pai1::keygen(PaiPublickey& paipk, PaiSecretkey& paisk)
{
    Pai::keygen(paipk, paisk);

	// Generate sk alpha from [1, lam]
    /*mpz_urandomm(alpha, rands, lam);
	mpz_mod(temp1, lam, alpha);
	while (mpz_cmp_ui(temp1, 0) != 0)
	{
		mpz_urandomm(alpha, rands, lam);
		mpz_mod(temp1, lam, alpha);		
	}*/
    
	mpz_divexact_ui (alpha, lam, 8);
	mpz_divexact(temp1, lam, alpha);
	mpz_powm(g, g, temp1, nn);

	
    // Pre-compute hp
    mpz_mul(temp2, p, p);
    mpz_powm(hp, g, alpha, temp2); // g^alpha % p^2
    mpz_sub_ui(hp, hp, 1);
    mpz_divexact(hp, hp, p);
    buffer = mpz_invert(hp, hp, p);
    
    
    // Pre-compute hq
    mpz_mul(temp3, q, q);
    mpz_powm(hq, g, alpha, temp3); // g^alpha % q^2
    mpz_sub_ui(hq, hq, 1);
    mpz_divexact(hq, hq, q);
    buffer = mpz_invert(hq, hq, q);
    
    
    // Precompute L(g^alpha % n^2)^-1
    mpz_powm(temp1, g, alpha, nn);
    mpz_sub_ui(temp1, temp1, 1);
    mpz_divexact(temp1, temp1, n);
    buffer = mpz_invert(precomp, temp1, n);
    
    
    // Set keys
    paipk.set (n, nn, g);
    paisk.set (p, q, temp2, temp3);
    paisk.set1(alpha);
}

inline PaiCiphertext& Pai1::encrypt(PaiPlaintext& paipt, PaiPublickey& paipk)
{
    // Get pt and pk
    paipt.get(m);
    paipk.get(n, nn, g);
    
    // Choose r from [0, n-1]
    mpz_urandomm(r, rands, n);
    
    // Compute ciphertext c = (g^(m+nr)) % n^2
    mpz_mul(temp1, r, n);
    mpz_add(temp1, temp1, m);
    mpz_powm(c, g, temp1, nn);
        
    // Set ciphertext
    ct.set(c);
    return ct;
}

inline PaiPlaintext& Pai1::decrypt(PaiCiphertext& paict, PaiPublickey& paipk, PaiSecretkey& paisk)
{
    // Get c, pk and sk
    paict.get(c);
    paisk.get1(alpha);
    paipk.get(n, nn);
    

    // Compute L(c^alpha % n^2)
    mpz_powm(temp1, c, alpha, nn);        
    mpz_sub_ui(temp1, temp1, 1);
    mpz_divexact(temp1, temp1, n);

    
    // Compute d = (L(.) * precomp) % n
    mpz_mul(temp1, temp1, precomp);
    mpz_mod(d, temp1, n);
    
    
    // Set plaintext
    pt.set(d);
    return pt;
}


inline PaiPlaintext& Pai1::decryptCRT(PaiCiphertext& paict, PaiPublickey& paipk, PaiSecretkey& paisk)
{
    // Get c, pk and sk
    paict.get(c);
    paisk.get1(alpha);
    paisk.get(p, q, temp2, temp3);
    
    
    // Compute mp
    mpz_powm(mp, c, alpha, temp2);
    mpz_sub_ui(mp, mp, 1);
    mpz_divexact(mp, mp, p);
    mpz_mul(mp, mp, hp);
    mpz_mod(mp, mp, p);
    
    
    // Compute mq
    mpz_powm(mq, c, alpha, temp3);
    mpz_sub_ui(mq, mq, 1);
    mpz_divexact(mq, mq, q);
    mpz_mul(mq, mq, hq);
    mpz_mod(mq, mq, q);
    
    
    // Compute m from mp and mq
    Pai::CRT(d, mp, mq, p, q);
    
    
    // Set plaintext
    pt.set(d);
    return pt;
}

#endif // Pai1_HPP
