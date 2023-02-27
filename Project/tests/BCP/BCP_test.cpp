/**
 * AUTHOR - Parthasarathi Das
 * BRIEF  - Measures the runtimes the BCP cryptosystem and writes timing information to files
 */


#include <iostream>
#include <fstream>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <gmp.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include "../../src/measurement.h"
#include "../../src/maxheaders.hpp"
#include "../../src/systems/BCP/BCP.hpp"
#include <vector>

NTL_CLIENT

int main(int argc, char** argv)
{
    int               mesbits         = strtol(argv[1], NULL, 10);
    int               mbits           = strtol(argv[2], NULL, 10);
    char              *type           = argv[3];
    
    int               ebits           = 0;
    int               parameter_count = 10;
    int               message_count   = 100;
    char              m_file[PATH_MAX], N_file[PATH_MAX], p_file[PATH_MAX], q_file[PATH_MAX], time_file[PATH_MAX];
    
    BCP b;
    
    ofstream          timeout;
    RR                factor;
    mpz_t             N, N2, m, c1, c2, dm, mod, twom, prime1, prime2, scalar, scalarm;
    gmp_randstate_t   rands;
    BCPPublickey      pk;
    BCPSecretkey      sk;
    BCPPlaintext      pt, *dt, *dtsum, *dtscal;
    BCPCiphertext     *ct, ct1, ct2, *ctsum, *ctscal;
    
    mpz_inits             (N, N2, m, c1, c2, dm, mod, twom, prime1, prime2, scalar, scalarm, NULL);
    gmp_randinit_default  (rands);
    
    
    
    // Open file for writing
    snprintf(time_file, PATH_MAX, "ED/%s_%s.txt", argv[1], argv[2]);
    timeout.open(time_file, std::ios_base::app);
    
    
    // Initialise variables
    if (strcmp(type, "short") == 0)
    {
        switch(mbits)
        {
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
                cout << "Error: Invalid security level\n";
                cout << "Usage: Valid inputs are 3072, 7680 and 15360\n";
                exit(0);
                break;
        }
    }
    else
    {
        ebits = 0;
    }
    
    
    // Set the name and path of .txt files for reading and  writing
    FILE *inm;
    FILE *inN;
    FILE *inp;
    FILE *inq;
    
    snprintf(m_file, PATH_MAX, "../../inputs/BP/m/%s.txt", argv[1]);
    snprintf(N_file, PATH_MAX, "../../inputs/BP/%s/%s_N.txt", argv[2],argv[2]);
    snprintf(p_file, PATH_MAX, "../../inputs/BP/%s/%s_p.txt", argv[2],argv[2]);
    snprintf(q_file, PATH_MAX, "../../inputs/BP/%s/%s_q.txt", argv[2],argv[2]);
    
    
    // DS for storing timing information
    Vec<ZZ> ticks;
    Vec<RR> runtime;
    ticks.SetLength(2);
    runtime.SetLength(2);
    
    
    // open message and N files
    inm  = fopen(m_file, "r");
    inN  = fopen(N_file, "r");
    inp  = fopen(p_file, "r");
    inq  = fopen(q_file, "r");
    
    for (int temp_count = 0; temp_count < parameter_count; temp_count++)
    {
        // Read N
        mpz_inp_str(mod, inN, 16);
        mpz_inp_str(prime1, inp, 16);
        mpz_inp_str(prime2, inq, 16);
        
        
        b.init(mbits, ebits, mod, prime1, prime2);
        b.keygen(pk, sk);
        
        // Encrypt and decrpyt for <message_count> messages
        for (int mcount = 0; mcount < message_count; mcount++)
        {
            // Set N, plaintext
            mpz_inp_str(m, inm, 16);
            pt.set(m);
            
            
            // Measure ED runtimes
            MEASURE ( ct = &b.encrypt(pt, pk); );
            ticks[0] = ticks[0] + (etime);
            
            MEASURE ( dt = &b.decrypt(*ct, pk, sk); );
            ticks[1] = ticks[1] + (etime);
            
            
            dt->get(dm);
            assert (mpz_cmp(m, dm) == 0);
        
        
        } // end-of-message-count-loop
        
    } // end-of-parameter-for-loop
    
        cout << mesbits << "\t" << mbits << "\t" << endl;
    
    // Compute factor before writing timing information to files
    factor = 1000 / (frequency * GHz * parameter_count * message_count);
    
    
    // Write timing information to file
    for (int k = 0; k < 2; k++)
    {
        runtime[k] = conv<RR>(ticks[k]) * factor;
        timeout << " & " << round(runtime[k]);
    }
    
    
    fclose(inm);
    fclose(inN);
    fclose(inp);
    fclose(inq);


    if (strcmp(type, "full") == 0)
        timeout << "\\\\" << endl;
    
    timeout.close();
    
    mpz_clears    (N, N2, m, c1, c2, dm, mod, prime1, prime2, twom, scalar, scalarm, NULL);
    gmp_randclear (rands);
    
    return 0;
}
