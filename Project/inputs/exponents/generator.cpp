/**
 
 @author Parthasarathi Das
 @desc   Generates messages of a given bit length
 
 */


#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <ntl/ZZ.h>
#include <ntl/ZZ_limbs.h>

NTL_CLIENT

void getB(mpz_t B, const mpz_t D, const int fbits)
{
    int i;
    mpfr_t  b, d, pi;
    
    mpfr_init2 (b,  200);
    mpfr_init2 (d,  200);
    mpfr_init2 (pi, 200);
    
    i = mpfr_const_pi  (pi,        MPFR_RNDU);    // Set pi value to pi
    i = mpfr_set_z     (d , D,     MPFR_RNDU);    // Get DK as mpfr_t type
    i = mpfr_log       (b , d,     MPFR_RNDU);    // 1. Compute log(|DK|)
    i = mpfr_sqrt      (d , d,     MPFR_RNDU);    // 2. Compute |DK|^1/2
    i = mpfr_mul       (b , b, d , MPFR_RNDU);    // 3. Multiply b & d
    i = mpfr_div       (b , b, pi, MPFR_RNDU);    // 4. Divide b by pi
    i = mpfr_div_ui    (b , b, 4 , MPFR_RNDU);    // 5. Divide b by 4
    i = mpfr_get_z     (B , b,     MPFR_RNDU);    // Set B = b
    
    mpz_mul_2exp(B, B, fbits); // Set B = 2^fbits * B
    
    mpfr_clear (b);
    mpfr_clear (d);
    mpfr_clear (pi);
    mpfr_free_cache ();
}

int main(int argc, char** argv)
{
    //
    int             fbits, Dbits, count, type;
    char            x_file[PATH_MAX], r_file[PATH_MAX];
    mpz_t           B, f, D;
    mpfr_t          b, d, pi;
    size_t          buffer;
    FILE            *outx, *outr;
    ZZ              x, r, UB;
    ofstream        xout, rout;
    gmp_randstate_t rands;
    
    // Initialise
    mpz_inits (B, f, D, NULL);
    gmp_randinit_default (rands);
    
    type   = strtol(argv[1], NULL, 10);
    if (argc != 2)
    {
        fbits  = strtol(argv[2], NULL, 10);
        Dbits  = strtol(argv[3], NULL, 10);
    }
    count  = 1000;
    
    
    // Set output file
    if (type == 0)
    {
        snprintf(x_file, PATH_MAX, "x_%d_%d.txt", fbits, Dbits);
        snprintf(r_file, PATH_MAX, "r_%d_%d.txt", fbits, Dbits);
        mpz_urandomb(f, rands, fbits);
        mpz_urandomb(D, rands, Dbits);
        getB(B, D, fbits);
        mpz_mul(B, B, f);
        ZZ_limbs_set(UB, B[0]._mp_d, B[0]._mp_size);
    }
    
    if ((type == 256) || (type == 384) || (type == 512))
    {
        snprintf(x_file, PATH_MAX, "x_%d.txt", type);
        snprintf(r_file, PATH_MAX, "r_%d.txt", type);
    }
    
    //
    xout.open(x_file);
    rout.open(r_file);
    
    // Generate x and r exponents
    if (type == 0)
    {
        for (int i = 0; i < count; i++)
        {
            x = RandomBnd(UB);
            r = RandomBnd(UB);
            xout << x << endl;
            rout << x << endl;
        }
    }
    
    else
    {
        for (int i = 0; i < count; i++)
        {
            x = RandomLen_ZZ(type);
            r = RandomLen_ZZ(type);
            xout << x << endl;
            rout << x << endl;
        }
    }
    
    
    // Close output files
    xout.close();
    rout.close();
    
    mpz_clears (B, f, D, NULL);
    gmp_randclear(rands);
    
    return 0;
}
