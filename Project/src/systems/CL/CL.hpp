/**
 *
 * AUTHOR - Parthasarathi Das
 * BRIEF - A class for the CL cryptosystem
 *
 */

#ifndef CL_HPP
#define CL_HPP

#include "../../parameters/CL/CLPublickey.hpp"
#include "../../parameters/CL/CLSecretkey.hpp"
#include "../../parameters/CL/CLPlaintext.hpp"
#include "../../parameters/CL/CLCiphertext.hpp"
#include "../../QFExponentiation.hpp"


NTL_CLIENT

class CL
{
    public:
    CL  ();
    CL  (const int window, const int mbits, const int dkbits, const int expbits);
    ~CL ();

    void initialise (const int no, const int exp, const mpz_t conductor, const mpz_t* arr, const mpz_t D);
    void keygen ();
    
    void solve (mpz_t dm, QF& M, const mpz_t f);
    void solve (mpz_t rop, QF& GG, const QF& GE, const mpz_t f, const mpz_t prime, const int exp); // PPH
    void solve (mpz_t rop, QF& F, QF& M, const mpz_t f, const mpz_t* pit, const mpz_t* pN, const int N, const int t);
    
    void CRT  (mpz_t x, mpz_t a, mpz_t b, mpz_t p1, mpz_t p2);
    void CRTl (mpz_t x, const int N, mpz_t *v, mpz_t *moduli);
    void getB (const mpz_t D, const int fbits);
    
    static int check_legendre (mpz_t* pN, const mpz_t p, const int i, const int symbol);
    
    
    
    
    
    
    protected:
    int                       prime_index;
    unsigned int              N, t, fbits, fmaxbits, DKbits, ebits, qbits, prime_r, temp;
    mpz_qform_t               form;
    mpz_qform_group_t         group;
    gmp_randstate_t           rands;
    mpz_t                     p, prod, q, r, k, B, UB, x, f, DK, Df, Df2, seed, m, m1, m2, dm, dm1, dm2;
    mpz_t                     *pN, *pit, *var, *mN;
    ZZ                        ZZ_p, ZZ_m, ZZ_x, ZZ_f, ZZ_UB;
    QF                        F, liftF, F1, F2, G, H, psiH, R, RR, psi_RR, C1, C2, C3, M, M1;
    QFExponentiation          qfe, qfe2, qfe3, qfed, qfed2, qfeF, qfeF1, qfeF2, qfeG, qfeH, qfeM, qfepsiH;
    vector<QF>                FN, lFN, TN, TN1, TN2, TN3, CN, MN;
    vector<QFExponentiation>  qfeFN, qfelFN, qfeTN;
    CLPlaintext               ptt;
    CLCiphertext              ctt;
    int len;
    
    // Buffer variables
    mpz_t  temp1, temp2, temp3, temp4, temp5, temp6;
    QF     T1, T2, T3, T4, T5, T6, CIX;
    ZZ     e1, e2, e3, zrop, te;
            
    // Variables for setB()
    mpfr_t b, d, pi;

    // Variables for solve()
    mpz_t st, stemp1, stemp2;
    
    // Variables for Chinese Remainder-ing
    mpz_t crta, crtb, crts, crtt, crtq, crtr;

}; // end of class CL

inline CL::CL()
{
}

inline CL::CL(const int window, const int mbits, const int dkbits, const int expbits)
{
    fbits    = mbits;
    ebits    = expbits;
    DKbits   = dkbits;
    fmaxbits = (dkbits/2) - 2;
    len      = 30;

    gmp_randinit_default  (rands);
    mpz_qform_init        (&group, &form);
    mpz_qform_group_init  (&group);
    
    mpz_inits (p, prod, q, r, k, B, UB, x, f, seed, DK, Df, Df2, m, m1, m2, dm, dm1, dm2, NULL);
    mpz_inits (temp1, temp2, temp3, temp4, temp5, temp6, NULL);
    mpz_inits (crta, crtb, crts, crtt, crtq, crtr, NULL);
    

    // Initialize mpz arrays
    pN   = mpz_init_array (len);
    mN   = mpz_init_array (len);
    pit  = mpz_init_array (len);
    var  = mpz_init_array (len);
    
    
    //
    FN     .resize (len);
    TN     .resize (len);
    TN1    .resize (len);
    TN2    .resize (len);
    TN3    .resize (len);
    CN     .resize (len);
    MN     .resize (len);
    lFN    .resize (len);
    qfeFN  .resize (len);
    qfeTN  .resize (len);
    qfelFN .resize (len);

    
    // Initialize Ideal arrays
    for (int i = 0; i < len; i++)
    {
        FN [i] .init (group);
        TN [i] .init (group);
        TN1[i] .init (group);
        TN2[i] .init (group);
        TN3[i] .init (group);
        CN [i] .init (group);
        MN [i] .init (group);
        lFN[i] .init (group);
    }
    

    //
    R      .init   (group);
    RR     .init   (group);
    psiH   .init   (group);
    psi_RR .init   (group);
    liftF  .init   (group);
    C1     .init   (group);
    C2     .init   (group);
    C3     .init   (group);
    CIX    .init   (group);
    M      .init   (group);
    M1     .init   (group);
    F      .init   (group);
    F1     .init   (group);
    F2     .init   (group);
    G      .init   (group);
    H      .init   (group);
    T1     .init   (group);
    T2     .init   (group);
    T3     .init   (group);
    T4     .init   (group);
    T5     .init   (group);
    T6     .init   (group);
    ctt    .init   (group);
    
    
    //    
    qfe   .init  (window);
    qfe2  .init  (window);
    qfe3  .init  (window);
    qfed  .init  (window);
    qfed2 .init  (window);
    qfeF  .init  (window);
    qfeF1 .init  (window);
    qfeF2 .init  (window);
    qfeG  .init  (window);
    qfeH  .init  (window);
    qfeM  .init  (window);
    
    qfepsiH.init(window);
    
    for (int i = 0; i < len; i++)
    {
        qfeFN [i] .init (window);
        qfeTN [i] .init (window);
        qfelFN[i] .init (window);
    }
    
    
    // setB temps
    mpfr_init2 (b,  200);
    mpfr_init2 (d,  200);
    mpfr_init2 (pi, 200);

    
    // solve4
    mpz_inits (st, stemp1, stemp2, NULL);
}



inline CL::~CL()
{
    gmp_randclear          (rands);
    mpz_qform_clear        (&group, &form);
    mpz_qform_group_clear  (&group);
    mpz_clears             (p, prod, q, r, k, B, UB, x, f, seed, DK, Df, Df2, m, m1, m2, dm, dm1, dm2, NULL);
    mpz_clears             (temp1, temp2, temp3, temp4, temp5, temp6, NULL);
    mpz_clears             (crta, crtb, crts, crtt, crtq, crtr, NULL);
    mpz_clear_array        (pN, len);
    mpz_clear_array        (mN, len);
    mpz_clear_array        (pit, len);
    mpz_clear_array        (var, len);
    

    mpfr_clear (b);
    mpfr_clear (d);
    mpfr_clear (pi);
    mpfr_free_cache ();

    mpz_clears (st, stemp1, stemp2, NULL);
}




inline void CL::initialise(const int no, const int exp, const mpz_t conductor, const mpz_t* arr, const mpz_t D)
{
    // Set no of primes in conductor and the power in the conductor
    N = no;
    t = exp;
    
    // Set conductor mpz_t and ZZ type
    mpz_set(f, conductor);
    ZZ_limbs_set(ZZ_f, f[0]._mp_d, f[0]._mp_size);
    
    // Set fundamental discriminant
    mpz_set(DK, D);
    
    // Set primes and prime powers of conductor
    for (int i = 0; i < N; i++)
    {
        mpz_set(pN[i], arr[i]);
        mpz_pow_ui(pit[i], pN[i], t);
    }

    // Set p if conductor is a single prime power
    if (N == 1)
        mpz_set(p, pN[0]);
}

inline void CL::keygen ()
{
    // Get upper bound on odd part; sets CL variable B
    getB(DK, fbits);
    
    // Set DK to -DK
    mpz_neg(DK, DK);
    
    // Compute Df
    mpz_mul(Df, f, f);
    mpz_mul(Df, Df, DK);
    
    // Compute Df2
    mpz_mul(Df2, f, f);
    mpz_mul(Df2, Df2, Df);
    
    // Set F to Castagnos ideal (f^2, f) of discriminant Df
    mpz_mul(temp1, f, f);
    F.assign(temp1, f, Df);
    
    // Compute ideal R in DK
    mpz_set_ui(r, 2);
    mpz_gcd(temp1, r, f);
    while ( !( (mpz_legendre(DK, r) == 1) && (mpz_cmp_ui(temp1, 1) == 0)))
    {
        mpz_nextprime(r, r);
        mpz_gcd(temp1, r, f);
    }
    prime_r = mpz_get_ui(r);
    prime_index = prime_index_ge(prime_r);
    prime_index = R.next_prime_QF(prime_index, DK);
    
    // Compute ideal R^2
    sqr(RR, R);
    
    // Generate a random k in (Z/conZ)^*
    mpz_urandomm (k, rands, f);
    mpz_gcd(temp1, k, f);
    while ((mpz_cmp_ui(k, 0) == 0) && (mpz_cmp_ui(temp1, 1) != 0))
    {
        mpz_urandomm (k, rands, f);
        mpz_gcd(temp1, k, f);
    }
    
    // Compute phi_inv(R^2)
    liftQF(psi_RR, RR, DK, Df, f);
    
    
    // PRECOMPUTATIONS
    
    e1 = RandomBnd(ebits);
    
    // Set Fi's for N > 1
    for (int i = 0; i < N; i++)
    {
        mpz_mul(temp1, pit[i], pit[i]);
        FN[i].assign(temp1, pit[i], Df);
        qfeFN[i].initialize(FN[i], e1);
    }
    
    
    // lift(F)
    liftQF(T2, F, Df, Df2, f);
    qfe.initialize(T2, ZZ_f);
    qfe.power(liftF, T2, ZZ_f);
    F.setD(Df);
    
    
    // Lift Fi's
    for (int i = 0; i < N; i++)
    {
        liftQF(T2, FN[i], Df, Df2, f);
        qfelFN[i].initialize(T2, ZZ_f);
        qfelFN[i].power(lFN[i], T2, ZZ_f);
        F.setD(Df);
    }
}

/**
 *
 * Computes an upper bound B on odd part of the class number s
 * Source - Page 12, https://eprint.iacr.org/2015/047.pdf
 *
 */
inline void CL::getB(const mpz_t D, const int fbits)
{
    int i;
    
    i = mpfr_const_pi(pi, MPFR_RNDU);       // Set pi value to pi
    i = mpfr_set_z(d, D, MPFR_RNDU);        // Get DK as mpfr_t type
    i = mpfr_log(b, d, MPFR_RNDU);          // 1. Compute log(|DK|)
    i = mpfr_sqrt(d, d, MPFR_RNDU);         // 2. Compute |DK|^1/2
    i = mpfr_mul(b, b, d, MPFR_RNDU);       // 3. Multiply b & d
    i = mpfr_div(b, b, pi, MPFR_RNDU);      // 4. Divide b by pi
    i = mpfr_div_ui(b, b, 4, MPFR_RNDU);    // 5. Divide b by 4
    i = mpfr_get_z(B, b, MPFR_RNDU);        // Set B = b
    
    mpz_mul_2exp(B, B, fbits);              // Set B = 2^fbits * B
}


/**
 *
 * Checks if the Legendre symbols (pN[j]/p) == (p/pN[j]) == leg
 * Returns 1 if true, 0 otherwise
 *
 */
inline int CL::check_legendre(mpz_t* pN, const mpz_t p, const int i, const int symbol)
{
    int j = 0;
    while (j < i)
    {
        if (!((mpz_legendre(pN[j], p) == symbol) && (mpz_legendre(p, pN[j]) == symbol)))
        {
            return 0;
        }

        j++;
    }
    return 1;
}


/**
 *
 * Solves the two simultaneous congruences, x = a (mod p1) and x = b (mod p2)
 * Source - http://doc.sagemath.org/html/en/reference/rings_standard/sage/arith/misc.html
 *
 */
inline void CL::CRT(mpz_t x, mpz_t a, mpz_t b, mpz_t p1, mpz_t p2)
{
    mpz_gcdext(x, crts, crtt, p1, p2);
    mpz_sub(crtt, b, a);
    mpz_fdiv_qr(crtq, crtr, crtt, x);
    mpz_lcm(crtt, p1, p2);
    mpz_mul(x, crtq, crts);
    mpz_mul(x, x, p1);
    mpz_add(x, x, a);
    mpz_mod(x, x, crtt);
}


/**
 *
 * Solves the simultaneous congruences, x = ai (mod pi)
 * Source - http://doc.sagemath.org/html/en/reference/rings_standard/sage/arith/misc.html
 *
 */
inline void CL::CRTl(mpz_t x, const int N, mpz_t *v, mpz_t *moduli)
{
    mpz_set (temp1, v[0]);
    mpz_set (temp2, moduli[0]);
    
    for (int i = 1; i < N; i++)
    {
        CRT (x, temp1, v[i], temp2, moduli[i]);
        mpz_lcm (temp2, temp2, moduli[i]);
        mpz_set (temp1, x);
    }
    mpz_mod(x, x, temp2);
}


/**
 *
 * Computes message m from reduced ideal M = (f^2, L(m)f)
 * Source - https://eprint.iacr.org/2015/047.pdf
 *
 */
inline void CL::solve (mpz_t dm, QF& M, const mpz_t f)
{
    // Get b part of M = (a,b)
    M.getb(dm);
        
    // Get m from b by computing multiplicative inverse modulo f
    mpz_divexact(dm, dm, f);
    mpz_invert(dm, dm, f);
}



inline void CL::solve(mpz_t rop, QF& GG, const QF& GE, const mpz_t f, const mpz_t prime, const int t)
{
    int buffer;
    
    // Initialise
    ZZ_limbs_set(ZZ_p, prime[0]._mp_d, prime[0]._mp_size);
    inv(T6, GG);
    qfe.initialize(GG, ZZ_p);
    qfe2.initialize(GE, ZZ_p);
    qfed2.initialize(T6, GE, ZZ_p, ZZ_p);
    
    // Initialise x_0 to 0
    mpz_set_ui(rop, 0);

	// For modulo p inversion
	mpz_set (temp2, prime);
    
    // Set variables for use when computing g^0
    GG.getD(temp1);
    mpz_set_ui(temp3, 1);
    
    // STEP 3:
    for(int k = 0; k < t; k++)
    {
        // Compute (g^(-x_k)*h)^(p^(t-1-k))
        buffer = t - 1 - k;
        e1 = power(ZZ_p, buffer);                          // Compute p^(e-1-k)
        ZZ_limbs_set(zrop, rop[0]._mp_d, rop[0]._mp_size); // Set (x_k)
        e2 = zrop * e1;                                    // Compute (x_k) (p^(e-1-k))
        
        if (zrop == 0)
            qfe2.power(T2, GE, e1);
        else
            qfed2.power(T2, T6, GE, e2, e1);
        
        T2.getb(temp3); // Extract b from ideal T2 = (a,b,c)
        
        // LOOP Step 2: Compute d_k(temp4)
        if (mpz_cmp_ui(temp3, 1) != 0)
        {
            mpz_divexact(temp4, temp3, temp2);
            mpz_invert(temp4, temp4, prime);
        }
        else
            mpz_set(temp4, prime);
        
        mpz_pow_ui(temp1, prime, k);  // Compute p^k
        mpz_mul(temp1, temp1, temp4); // Compute (p^k)(d_k)
        
        // LOOP STEP 3: Compute x_k + 1
        mpz_add(rop, rop, temp1);
    }
    mpz_mod(rop, rop, f);
}





/**
 *
 * Computes message m from ideal M when f = (p1p2...pN)^t
 * Uses generic Pohlig Hellman subroutine
 * Pseudocode/Algorithm source - https://en.wikipedia.org/wiki/Pohlig%E2%80%93Hellman_algorithm
 *
 */
inline void CL::solve(mpz_t rop, QF& F, QF& M, const mpz_t f, const mpz_t* pit, const mpz_t* pN, const int N, const int t)
{
    int temp;
    qfeM.initialize(M, e1);
    
    //Initial value
    mpz_set_ui(rop, 0);

    //
    for(int i = 0; i < N; i++)
    {
        mpz_divexact(stemp2, f, pit[i]);
        ZZ_limbs_set(e1, stemp2[0]._mp_d, stemp2[0]._mp_size);
        
        qfeF.power(T3, F, e1);
        qfeM.power(T4, M, e1);
        solve(st, T3, T4, pit[i], pN[i], t);
        mpz_set(var[i], st);
        
        // CRT steps
        temp = mpz_invert(st, stemp2, pit[i]);
        mpz_mul(st, st, stemp2);
        mpz_mul(st, st, var[i]);
        mpz_add(rop, rop, st);
    }
    
    mpz_mod(rop, rop, f);
}

#endif // CL_HPP
