/**
 * AUTHOR - Parthasarathi Das
 * BRIEF  - Measures the runtimes the Paillier cryptosystem (basic) and writes timing information to files
 */

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <string.h>
#include <unistd.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include "../../../src/measurement.h"
#include "../../../src/maxheaders.hpp"
#include "../../../src/systems/Pai/Pai.hpp"
#include "../../../src/systems/Pai/Pai0.hpp"
#include <vector>

NTL_CLIENT

int main(int argc, char** argv)
{
    int               mesbits         = strtol(argv[1], NULL, 10);
    int               mbits           = strtol(argv[2], NULL, 10);
    
    int               parameter_count = 1;
    int               message_count   = 100;
    
    char              m_file[PATH_MAX];
    char              N_file[PATH_MAX];
    char              p_file[PATH_MAX];
    char              q_file[PATH_MAX];
    char              t_file[PATH_MAX];
    
    Pai0 b;
    
    ofstream          timeout;
    RR                factor;
    mpz_t             N, N2, p, q, m, dm, c1, c2, mod, twom, scalar, scalarm;
    gmp_randstate_t   rands;
    PaiPublickey      pk;
    PaiSecretkey      sk;
    PaiPlaintext      pt, *dt, *dtsum, *dtscal;
    PaiCiphertext     *ct, ct1, ct2, *ctsum, *ctscal;
    
    mpz_inits            (N, N2, p, q, m, dm, c1, c2, mod, twom, scalar, scalarm, NULL);
    gmp_randinit_default (rands);
    
    
    // Set the name and path of .txt files for reading and  writing
    FILE *inm;
    FILE *inN;
    FILE *inp;
    FILE *inq;

    snprintf(m_file, PATH_MAX, "../../../inputs/BP/m/%s.txt", argv[1]);
    snprintf(N_file, PATH_MAX, "../../../inputs/BP/%s/%s_N.txt", argv[2], argv[2]);
    snprintf(p_file, PATH_MAX, "../../../inputs/BP/%s/%s_p.txt", argv[2], argv[2]);
    snprintf(q_file, PATH_MAX, "../../../inputs/BP/%s/%s_q.txt", argv[2], argv[2]);
    snprintf(t_file, PATH_MAX, "ED/%s_%s.txt", argv[1], argv[2]);

    
    // Open file for writing
    timeout.open(t_file);
    
    // DS for storing timing information
    Vec<ZZ> ticks;
    Vec<RR> runtime;
    ticks.SetLength(3);
    runtime.SetLength(3);
    
    
    // open message file
    inm  = fopen(m_file,  "r");
    inN  = fopen(N_file,  "r");
    inp  = fopen(p_file,  "r");
    inq  = fopen(q_file,  "r");
    
    // Start
    for (int temp_count = 0; temp_count < parameter_count; temp_count++)
    {
        mpz_inp_str(mod, inN, 10);
        mpz_inp_str(p,   inp, 10);
        mpz_inp_str(q,   inq, 10);
        
        b.init(mbits, mod, p, q);
        b.keygen(pk, sk);
        
        // Encrypt and decrpyt for message_count messages
        for (int mcount = 0; mcount < message_count; mcount++)
        {
            // Set plaintext
            mpz_inp_str(m, inm, 16);
            pt.set(m);
            
            
            // Measure function runtimes
            MEASURE ( ct = &b.encrypt(pt, pk); );
            ticks[0] = ticks[0] + (etime);
            
            MEASURE ( dt = &b.decrypt(*ct, pk, sk); );
            ticks[1] = ticks[1] + (etime);
            
            dt->get(dm);
            assert (mpz_cmp(m, dm) == 0);
            mpz_set_ui(dm, 0);
            
            MEASURE ( dt = &b.decryptCRT(*ct, pk, sk); );
            ticks[2] = ticks[2] + (etime);
            
            dt->get(dm);
            assert (mpz_cmp(m, dm) == 0);
            
        } // end-of-message-count-loop
        
    } // end-of-parameter-for-loop
    
    cout << mesbits << " " << mbits << endl;
    
    // Compute factor before writing timing information to files
    factor = 1000 / (frequency * GHz * parameter_count * message_count);
    
    
    // Write timing information to file
    for (int k = 0; k < 3; k++)
    {
        runtime[k] = conv<RR>(ticks[k]) * factor;
        timeout << " & " << round(runtime[k]);
    }
    
    
    timeout.close();
    fclose(inm);
    fclose(inN);
    fclose(inp);
    fclose(inq);
    
    mpz_clears (N, N2, p, q, m, dm, c1, c2, mod, twom, scalar, scalarm, NULL);
    gmp_randclear (rands);
    
    return 0;
}
