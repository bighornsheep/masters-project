#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <mpfr.h>
#include <NTL/ZZ.h>
#include <iostream>
#include <fstream>
#include <vector>

#include "../../src/measurement.h"
#include "../../src/maxheaders.hpp"
#include "../../src/QF.hpp"
#include "../../src/systems/CL/Variant.hpp"

NTL_CLIENT


// Forward declarations
void txt_write (const int N, const int t, const int fbits, const Vec<RR> edtime, std::ofstream &redout);
void csv_write (const int N, const int t, const int fbits, const Vec<RR> edtime, std::ofstream &redout);


int main(int argc, char** argv)
{
    int               dkbits          = strtol(argv[1], NULL, 10);
    int               fbits           = strtol(argv[2], NULL, 10);
    char              *type           =        argv[3];
    int               parameter_count = strtol(argv[4], NULL, 10);
    int               message_count   = strtol(argv[5], NULL, 10);

    int               ebits           =  0;
    int               fmaxbits        =  0;
    int               count           =  0;
    int               mcounter        =  0;
    int               window_size     =  8;
    int               N               = 24;
    int               t               =  5;

    
    char              infile_m  [PATH_MAX];
    char              infile_f  [PATH_MAX];
    char              infile_fa [PATH_MAX];
    char              infile_pr [PATH_MAX];
    char              infile_dk [PATH_MAX];
    
    CLPublickey       clpk;
    CLSecretkey       clsk;
    CLPlaintext       clpt, *cldt, *cldtsum, *cldtscal;
    CLCiphertext      clct1, clct2, *clct, *clctsum, *clctscal;
    
    mpz_t             m, dm, twom, alpha, alpham, c1, c2, c3, c4, DK, prod, con, gcd, seed, temp, *pN;
    RR                factor;
    QF                CC1, CC2;
    vector<QF>        CC1N, CC2N;
    mpz_qform_t       form1;
    gmp_randstate_t   rands;
    mpz_qform_group_t group1;
    
    
    // Initialise variables
    if (strcmp(type, "short") == 0)
    {
        switch(dkbits)
        {
            case(1828):
                ebits = 2 * 128;
                break;
            case(3598):
                ebits = 2 * 192;
                break;
            case(5972):
                ebits = 2 * 256;
                break;

            case(3072):
                ebits = 2 * 128;
                break;
            case(7680):
                ebits = 2 * 192;
                break;
            case(15360):
                ebits = 2 * 256;
                break;

            default:
                cout << "Error: Check inputs\n";
                exit(0);
                break;
        }
    }

    else if (strcmp(type, "full") == 0) {}
    
    else
    {
        cout << "Error: Incorrect exponent size argument\n";
        exit(0);
    }
    
    //
    mpz_inits (m, dm, DK, prod, con, gcd, seed, temp, NULL);
    mpz_inits (c1, c2, c3, c4, twom, alpha, alpham, NULL);
    mpz_qform_init (&group1, &form1);
    mpz_qform_group_init (&group1);
    gmp_randinit_default (rands);
    pN = mpz_init_array(N);
    CC1.init (group1);
    CC2.init (group1);
    clct1.init (group1);
    clct2.init (group1);
    CC1N.resize(N);
    CC2N.resize(N);
    for (int i = 0; i < N; i++)
    {
        CC1N[i].init (group1);
        CC2N[i].init (group1);
    }
    


    Variant obj(window_size, fbits, dkbits, ebits);



	// Open output file
	char outfile_red  [PATH_MAX];
	ofstream redout;
	snprintf (outfile_red, PATH_MAX, "tables/V_%s_%s_%s.txt", argv[1], argv[2], argv[3]);
    redout.open (outfile_red);


    
    // Set fmaxbits
    fmaxbits = (dkbits / 2) - 2;
    
    
    //
    for (int i = 1; i < N + 1; i++)
    {
        for (int j = 1; j < t + 1; j++)
        {
            snprintf (infile_m,  PATH_MAX, "../../inputs/messages/%s.txt",       argv[2]               );
            snprintf (infile_f,  PATH_MAX, "../../inputs/CL/%s/%d_%d_f.txt" ,    argv[2],          i, j);
            snprintf (infile_fa, PATH_MAX, "../../inputs/CL/%s/%d_%d_fa.txt",    argv[2],          i, j);
            snprintf (infile_pr, PATH_MAX, "../../inputs/CL/%s/%d_%d_pr.txt",    argv[2],          i, j);
            snprintf (infile_dk, PATH_MAX, "../../inputs/CL/%s/%s/%d_%d_dk.txt", argv[2], argv[1], i, j);
            
            
            Vec<ZZ> edticks;
            Vec<ZZ> scticks;
            Vec<RR> edtime;
            Vec<RR> sctime;
            
            edticks .SetLength (7);
            edtime  .SetLength (7);
            scticks .SetLength (6);
            sctime  .SetLength (6);
            

            // if file exists, begin timing computation
            if ( access (infile_f, F_OK) != -1 )
            {
                count = 0;
                for (int temp_count = 0; temp_count < parameter_count; temp_count++)
                {
                    FILE *inf  = fopen(infile_f,  "r");
                    FILE *infa = fopen(infile_fa, "r");
                    FILE *inpr = fopen(infile_pr, "r");
                    FILE *indk = fopen(infile_dk, "r");
                    
                    
                    // read inputs
                    mpz_inp_str(prod, inpr, 16);
                    mpz_inp_str(con, inf, 16);
                    mpz_inp_str(DK, indk, 16);
                    for (int k = 0; k < i; k++)
                        mpz_inp_str(pN[k], infa, 16);

                    
                    // send parameters to keygen
                    obj.initialise(i, j, con, pN, DK);
                    obj.keygen(clpk, clsk);
                    
                    
                    // encrypt and decrpyt each message
                    mcounter = 0;
                    for (int mcount = 0; mcount < message_count; mcount++)
                    {
                        FILE *inm  = fopen(infile_m,  "r");
                        
                        // Get message
                        mpz_inp_str(m, inm, 16);
                        
                        // Check gcd
                        mpz_gcd(gcd, m, con);
                        if ((mpz_cmp_ui(gcd, 1) == 0) && mpz_cmp(m, con) < 0)
                        {
                            mcounter++;
                            clpt.set(m);
                            
							if (i < 10)
							{
                            MEASURE ( clct = &obj.encrypt(clpt, clpk); );
                            edticks[0] = edticks[0] + (etime);
                            }

                            if (j == 1)
                            {
                                
                                // 1. Modular inverse decryption
                                MEASURE ( cldt = &obj.decrypt(*clct, clpk, clsk); );
                                edticks[1] = edticks[1] + (etime);
                                cldt->get(dm);
                                assert (mpz_cmp(m, dm) == 0);
                                
                                if (i > 1)
                                {
                                    // 1. CRT 1 decryption
                                    MEASURE ( cldt = &obj.decryptCRT(*clct, clpk, clsk); );
                                    edticks[2] = edticks[2] + (etime);
                                    cldt->get(dm);
                                    assert (mpz_cmp(m, dm) == 0);                                    
                                    
                                    
                                    // 2. CRT 2 encryption and decryption
                                    MEASURE ( clct = &obj.encrypt2(clpt, clpk); );
                                    edticks[3] = edticks[3] + (etime);
                                    MEASURE ( cldt = &obj.decrypt2(*clct, clpk, clsk); );
                                    edticks[4] = edticks[4] + (etime);
                                    cldt->get(dm);
                                    assert (mpz_cmp(m, dm) == 0);
                                    
                                    
                                    // 3. CRT 3 encryption and decryption
                                    MEASURE ( clct = &obj.encrypt3(clpt, clpk); );
                                    edticks[5] = edticks[5] + (etime);
                                    MEASURE ( cldt = &obj.decrypt3(*clct, clpk, clsk); );
                                    edticks[6] = edticks[6] + (etime);
                                    cldt->get(dm);
                                    assert (mpz_cmp(m, dm) == 0);
                                }
                            }
                            
                            
                            
                            else
                            {
								if (i < 10)
								{
                                // 1. Pohlig Hellman decryption
                                MEASURE ( cldt = &obj.hdecrypt(*clct, clpk, clsk); );
                                edticks[1] = edticks[1] + (etime);
                                cldt->get(dm);
                                assert (mpz_cmp(m, dm) == 0);
                                }
                                if (i > 1)
                                {
                                    // ENCRYPTION FOR CRT 3
                                    MEASURE ( clct = &obj.encrypt3(clpt, clpk); );
                                    edticks[5] = edticks[5] + (etime);
                                    
                                    // 1. CRT 3 in PH Decryption
                                    MEASURE ( cldt = &obj.hdecrypt3(*clct, clpk, clsk); );
                                    edticks[6] = edticks[6] + (etime);
                                    cldt->get(dm);
                                    assert (mpz_cmp(m, dm) == 0);
                                }
                            }
                            
                        } // gcd and m size check
                        
                        fclose(inm);
                        
                    } // end of message loop
                    
                    count = count + mcounter;
                    
                    // Close parameter files
                    fclose(inf);
                    fclose(infa);
                    fclose(inpr);
                    fclose(indk);
                    
                                        
                } // end of parameter loop
                
                // Print timing information on screen
                cout << count << "\t iterations for " << i << " " << j << endl;

                
                // Compute run time in milliseconds
                factor = 1000 / (frequency * GHz * count);
                for (int k = 0; k < 7; k++)
                {
                    edtime[k] = conv<RR>(edticks[k]) * factor;
                }

				// Write to text file
				txt_write (i, j, fbits, edtime, redout );
                
            } // end of file check condition

        } //j
    
    }//i

    // Close time files
    redout .close();

    
    // Deallocate memory
    mpz_clears            (m, dm, DK, prod, con, seed, gcd, temp, NULL);
    mpz_clears            (c1, c2, c3, c4, twom, alpha, alpham, NULL);
    mpz_clear_array       (pN, N);
    gmp_randclear         (rands);
    mpz_qform_clear       (&group1, &form1);
    mpz_qform_group_clear (&group1);
    
    return 0;
}




void txt_write (const int N, const int t, const int fbits, const Vec<RR> edtime, std::ofstream &redout)
{

	// Write raw E D timing information to file
    redout << N << " & " << t << " & " ;

	if ((fbits == 3072) or (fbits == 7680) or (fbits == 15360))
	{
		redout << round(edtime[0]) << " & ";
		redout << round(edtime[1]) << " & ";

		if (edtime[5] != 0) 
		{
			redout << round(edtime[5]) << " & ";
			redout << round(edtime[6]);
		}

		redout << " \\\\ " << endl;
	}

	else
	{
		for (int k = 0; k < 7; k++)
		{
			if (edtime[k] != 0)
			{
				if (k != 6)
					redout << round(edtime[k]) << " & ";
                else
                    redout << round(edtime[k]);
            }
                        
            else
            {
				if (k != 6)
                    redout << " & ";
            }
        }
        redout << " \\\\ " << endl;
    }
}

