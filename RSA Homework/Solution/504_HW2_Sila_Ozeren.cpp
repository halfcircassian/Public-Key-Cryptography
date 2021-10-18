/* 
;==========================================
; Title:  Implementation of RSA Algorithm 
; Author: Sila Ozeren
; Student ID: 2158087
; Date:   29 April 2021
; Operating System: Ubuntu
; IDE: Visual Studio Code
;==========================================
*/

#include <iostream>
#include <gmpxx.h>
#include <time.h>
#include <fstream>
#include <string>

using namespace std;

int main(){

    int pressChoice;

    std::cout << "******************MENU******************" << std::endl;
    std::cout << "To generate the RSA parameters, PRESS 1." << std::endl;
    std::cout << "To encrypt, PRESS 2." << std::endl;
    std::cout << "To decrypt, PRESS 3." << std::endl;
    std::cout << "******************MENU******************" << std::endl;

    std::cin >> pressChoice;

    if (pressChoice == 1)
    {
        gmp_printf("Processing...\n\n");
        
        gmp_randstate_t grt; 

        gmp_randinit_default (grt); 

        gmp_randseed_ui (grt, time (NULL)); 

        mpz_class key_p, key_q; 

    	mpz_urandomb (key_p.get_mpz_t(), grt, 512); 
    
        /*It is important to note that key_p is ranging over 0 - 2^512. I want to make
        sure that these prime parameters are exactly 512 bits so that 128 byte message 
        can be uniquely encrypted in mod n. Hence, */

	    mpz_class base = 2;
	    mpz_class exponent = 511;
        mpz_class OR;
	    mpz_pow_ui(OR.get_mpz_t(), base.get_mpz_t(), exponent.get_ui());
	    key_p = OR | key_p;

        mpz_nextprime (key_p.get_mpz_t(), key_p.get_mpz_t()); 

        mpz_urandomb(key_q.get_mpz_t(), grt, 512); // generates a random number a ^ 0-2 512
        mpz_pow_ui(OR.get_mpz_t(), base.get_mpz_t(), exponent.get_ui());
	    key_q = OR | key_q;
	
        mpz_nextprime (key_q.get_mpz_t(), key_q.get_mpz_t()); 

        gmp_printf("Prime p is:\n");
        gmp_printf("%ZX\n", key_p);


        int b = mpz_probab_prime_p(key_p.get_mpz_t(), 15); // Test if they are actually prime integers.

        //Return 2 if it is definitely prime, 
        //Return 1 if it is probably prime (without being certain), 
        //Return 0 if it is definitely non-prime.

        if (b == 2){
            gmp_printf("p is definitely a prime number.\n\n");
        }
        else if (b == 1){
            gmp_printf("p is probably a prime number.\n\n");
        }
        else if (b == 0){
            gmp_printf("p is definitely not a prime number.\n\n");
        }

        gmp_printf("Prime q is:\n");
	    gmp_printf("%ZX\n", key_q);

        //Return 2 if it is definitely prime, 
        //Return 1 if it is probably prime (without being certain),
        //Return 0 if it is definitely non-prime.
        
        int a = mpz_probab_prime_p(key_q.get_mpz_t(), 15); // reps = 15

        if (a == 2){
            gmp_printf("q is definitely a prime number.\n\n");
        }
        else if (a == 1){
            gmp_printf("q is probably a prime number.\n\n");
        }
        else if (a == 0){
            gmp_printf("q is definitely not a prime number.\n\n");
        }

        // generate the parameter n
        mpz_class key_n = key_p * key_q;

        // generate the phi of n
        mpz_class key_pMinus1 = key_p - 1;
        mpz_class key_qMinus1 = key_q - 1;

        mpz_class key_phiOfn = key_pMinus1 * key_qMinus1;

        gmp_printf("Phi(n) is:\n");
        gmp_printf("%ZX\n\n", key_phiOfn);

        // choose the public exponent e such that gcd(e,phiOfn) = 1
        // 65537 =  2^16 + 1 = (1 0000 0000 0000 0001)_{2}. #MUL + #SQ = 17. Hamming Weight is small.

        mpz_class key_e;
        mpz_init_set_ui(key_e.get_mpz_t(), 65537); 
        gmp_printf("Public exponent e is: \n");
        gmp_printf("%ZX\n\n", key_e);

        // generate the private key d such that e*d = 1(mod phiOfn)
        mpz_class key_d;
        mpz_invert(key_d.get_mpz_t(), key_e.get_mpz_t(), key_phiOfn.get_mpz_t());
        gmp_printf("Private key d is: \n");
        gmp_printf("%ZX\n\n", key_d);

        // saving to the file:
        ofstream parameters ("parameters.txt");
	    parameters << key_p << endl << key_q << endl << key_phiOfn <<endl << key_n << endl << key_e <<endl << key_d <<endl;
	    parameters.close();
    }

    else if (pressChoice == 2)
    {
        
        gmp_printf("Processing...\n\n");

        // reading required parameters. 
        mpz_class p, q, phiOfn, n, e, d;

        mpz_class numericCiphertext;

        ifstream parameters ("parameters.txt");
        while (parameters)
        parameters >> p >> q >> phiOfn >> n >> e >> d;
        parameters.close();

        /* "I say let the world go to hell, but I should always have my tea. :):):):):):):):)
        :):):):):):):):):):):):):):):):):):):):):):):)" is my message without the quotation marks. */

        FILE *file;
	    file = fopen("plain.txt", "r");

	    char arr[128];
	    fread(arr, sizeof(arr), 1, file);

        fclose (file);

	    mpz_class numericPlaintext;
	    size_t size = 128;
	    mpz_import(numericPlaintext.get_mpz_t(), 128, 1, sizeof(arr[0]), 0, 0, arr);

        cout << "Plaintext as a numeric value: " << numericPlaintext << endl << endl;

        // encryption;
	    mpz_powm(numericCiphertext.get_mpz_t(), numericPlaintext.get_mpz_t(), e.get_mpz_t(), n.get_mpz_t());

        cout << "Ciphertext as a numeric value: " << numericCiphertext << endl << endl;
        
	    mpz_export(arr, &size, 1, sizeof(arr[0]), 0, 0, numericCiphertext.get_mpz_t());

	    FILE *fp = fopen("cipher.txt", "w");
	    fwrite(arr, 1, 128, fp);

        fclose (fp);

    }

    else if (pressChoice == 3)
    {

        gmp_printf("Processing...\n\n");

        // reading required parameters. 
        mpz_class p, q, phiOfn, n, e, d, ct, pt;

        ifstream parameters ("parameters.txt");
        while (parameters)
        parameters >> p >> q >> phiOfn >> n >> e >> d;
        parameters.close();

        FILE *file;
	    file = fopen("cipher.txt", "r");

	    char arr[128];
	    fread(arr, sizeof(arr), 1, file);

        fclose (file);

	    mpz_class numCipher;
        mpz_class numPlain;

	    size_t size = 128;
	    mpz_import(numCipher.get_mpz_t(), 128, 1, sizeof(arr[0]), 0, 0, arr);

        // decryption
	    mpz_powm(numPlain.get_mpz_t(), numCipher.get_mpz_t(), d.get_mpz_t(), n.get_mpz_t());
	    mpz_export(arr, &size, 1, sizeof(arr[0]), 0, 0, numPlain.get_mpz_t());

	    FILE *fp = fopen("encodedMessage.txt", "w");
	    fwrite(arr, 1, 128, fp);

        fclose (fp);

        cout << "Decoded message is: " << arr;


    }  
    
    return 0;

    }

