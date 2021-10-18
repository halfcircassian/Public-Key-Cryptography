/* 
;==========================================
; Title:  Implementation of DHKE and Koblitz Imbedding 
; Author: Sila Ozeren
; Student ID: 2158087
; Date:   04 July 2021
; Operating System: Ubuntu
; IDE: Visual Studio Code
;==========================================
*/

/*
;==========================================
; I can ensure you that I put my heart and soul into this assignment. Since I had zero background in coding, it felt like a torture. 
; Nevertheless, I somehow mananged to finish it. The only problem is I could not read these parameters without manipulating the text
; file itself. Even though I could read a,b and p, I could not jump a blank line ect. So, I decided not to do so and take these 
; parameters by myself. Like a cave-woman. I know. I feel geniunly sorry for that. Also, I had to read the strings from three different
; files. Again, I am also sorry for that. If you consider a "file opening for dummies in c++" course, please sign me up. Thank you for this
; semestre. I have learnt a lot. Below, you can see my ec points.
;==========================================
*/

/*
;----------------------------------------------------------------------------------------------------------------------------------
;Is it imbed or embed?
;x coordinate of the message as an elliptic curve is: 824424009783958023608561221290166356779307564435202860
;y coordinate of the message as an elliptic curve is: 114603175179670441643583957226880527684408001728902527802917708171717045034557
;----------------------------------------------------------------------------------------------------------------------------------- 
;Koblitz calls it imbed.
;x coordinate of the message as an elliptic curve is: 55489894843533584385430682014313539037029998022270436336940
;y coordinate of the message as an elliptic curve is: 61838088649899776705462049437448151745188488494087866009825857113334982985832
;---------------------------------------------------------------------------------------------------------------------------------- 
;I need a vacation.
;x coordinate of the message as an elliptic curve is: 48923386754441804512954720853998389368126530863
;y coordinate of the message as an elliptic curve is: 104224915224219896085049881846866141044714890914718497374712262326506287072146
;----------------------------------------------------------------------------------------------------------------------------------
*/


#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <time.h>
#include <gmp.h>
#include <string.h>
#include <fstream>

using namespace std;

#define PRINT_TEXT_BLOCK 1
#if PRINT_TEXT_BLOCK
#define PrintFormattedTextBlock PrintFormattedTextBlock_
#else
#define PrintFormattedTextBlock(...)
#endif

void
PrintFormattedTextBlock_(const char *Prefix, uint32_t Size, char *Text)
{
  printf("character count: %" PRIu32 "\n", Size);
  printf("== %s text begin ==\n", Prefix);
  // printf will stop if you get a null in the string, use fwrite for dumping byte strings.
  fwrite(Text, 1, Size, stdout);
  printf("\n== %s text end ==\n", Prefix);
}

static void
EncodeFileIntoNumber(const char *Filename, mpz_t Out)
{
  FILE *File = fopen(Filename, "rb");
  if(File)
  {
    char *Buffer = 0;
    uint32_t FileSize;
    
    fseek(File, 0, SEEK_END);
    FileSize = ftell(File);
    rewind(File);
    
    Buffer = (char *)malloc((FileSize)*sizeof(char));
    if(Buffer)
    {
      fread(Buffer, 1, FileSize, File);
      PrintFormattedTextBlock("encoded", FileSize, Buffer);
      
      mpz_import(Out, FileSize, 1, 1, 0, 0, Buffer);
      free(Buffer);
    }
    fclose(File);
  }
}

static void
DecodeNumberIntoFile(const char *Filename, mpz_t In)
{
  FILE *File = fopen(Filename, "wb");
  if(File)
  {
    size_t StringSize = 0;
    
    char *Decoded = (char*)mpz_export(0, &StringSize, 1, 1, 0, 0, In);
    PrintFormattedTextBlock("decoded", (uint32_t)StringSize, Decoded);
    
    fwrite(Decoded, 1, StringSize, File);
    free(Decoded);
    fclose(File);
  }
}

void 
GF_addition(mpz_t result, mpz_t a, mpz_t b, mpz_t prime_p){

      mpz_t temp;
      mpz_init(temp);

      mpz_add(temp,a,b);

      /* Since we are dealing with ECC in a finite prime field, the result of a group operation
      must remain in the group, GF(p). 
  
      Check if temp is bigger than p. Comparision function mpz_cmp is used.*/

      if (mpz_cmp(temp, prime_p) >= 0){
          mpz_sub(result, temp, prime_p);

      /* 0UL is a c++ constant, which means unsigned long.*/

      }
      else if(mpz_cmp_ui(temp, 0UL) < 0){
          mpz_add(result, temp, prime_p);
      }
      else
          mpz_set(result, temp);

      mpz_clear(temp);
} 

/* This function is to perform substruction in a prime field GF(p). */ 
void 
GF_substruction(mpz_t result, mpz_t a, mpz_t b, mpz_t prime_p){

      mpz_t temp;
      mpz_init(temp);

      /* Negative of b is assigned to the variable mpz_t temp. */
      mpz_neg(temp, b); 

      GF_addition(result, a, temp, prime_p); 

      mpz_clear(temp);
}
    
struct 
Point{

      mpz_t a;
      mpz_t b;

};

/* This function is to create a new point at (0,0) */
struct 
Point *create_point(void)
{
    /*
    1) struct Point *new_point
    This statement declares a variable named "new_point" as type of pointer to struct Point. The "Point structure" was
    defined before this code snippet.

    2) malloc(sizeof(struct new_point));
    malloc() is a standard library function in C that allocates memory of size equal to the size of 
    new_point structure in the heap region during runtime and returns a void pointer to that location.
    A void pointer is a pointer that has no associated data type with it. A void pointer can hold address 
    of any type and can be typcasted to any type.
    
    3) new_point = (struct Point*) malloc(...); 
    The pointer returned by malloc() function is type-casted into pointer to struct Point type 
    and then assigned to new_point variable.
    
     */

	  struct Point *new_point = (struct Point*)malloc(sizeof(*new_point));
    
    /*To access members of a structure using pointers, we use the -> operator.*/

	  mpz_init_set_ui(new_point->a, 0UL); /*0UL = unsigned long with value 0*/
	  mpz_init_set_ui(new_point->b, 0UL);
	  return new_point; 
    
};

struct 
Point 
*copy_point(struct Point *point)
{
	struct Point *copy = create_point();
	mpz_set(copy->a, point->a);
	mpz_set(copy->b, point->b);
	return copy;
};

void 
clear_point(struct Point *point)
{
	mpz_clear(point->a);
	mpz_clear(point->b);
	free(point);
};

struct 
Elliptic_Curve{
    
    /*y^2 = x^3 + ax + b mod p*/

      mpz_t prime; 
      mpz_t a;     
      mpz_t b;     
      mpz_t n; 

  
};

/* This function if to perform point addition, where P != Q. 

    s = (Q_y - P_y) . (Q_x - P_x)^{-1} mod p; if P != Q (point addition)

    x_3 = s^2 - (P_x + Q_x) mod p, 
    y_3 = s(P_x - x_3) - P_y mod p

    For furher information, 
    Christof Paar · Jan Pelzl
    Understanding Cryptography
    A Textbook for Students and Practitioners

*/
struct 
Point *point_addition(struct Point *p, struct Point *q, struct Elliptic_Curve *ecc){


    /*If one of the points is (0,0), returns the non-zero point.*/
    if (mpz_cmp_ui(p->a, 0UL) == 0 && mpz_cmp_ui(p->b, 0UL) == 0)
		return copy_point(q);
	  if (mpz_cmp_ui(q->a, 0UL) == 0 && mpz_cmp_ui(q->b, 0UL) == 0)
		return copy_point(p);
  
    /*Create a point r at (0,0)*/
    struct Point *r = create_point(); 

    /*Check if P_x = Q_x = 0. If this is the case, slope would be undefined.*/
    if (mpz_cmp(p->a, q->a) == 0) 
		return r;

    /* Find Q_x - P_x mod p. */
    mpz_t x_alpha;
    mpz_init(x_alpha);
    GF_substruction(x_alpha, q->a, p->a, ecc->prime);

    /* Find Q_y - P_y mod p. */
    mpz_t y_alpha;
    mpz_init(y_alpha);
    GF_substruction(y_alpha, q->b, p->b, ecc->prime);

    /* Find (Q_x - P_x)^{-1} mod p. */

    mpz_t inverse_x_alpha;
    mpz_init(inverse_x_alpha);
    mpz_invert(inverse_x_alpha, x_alpha, ecc->prime);
   
    /* Perform  s = (Q_y - P_y) . (Q_x - P_x)^{-1} mod p */

    mpz_t s;
    mpz_init(s);
    mpz_mul(s, inverse_x_alpha, y_alpha);
    mpz_mod(s, s, ecc->prime);

    /*It is time to find x_3 and y_3.
    
    x_3 = s^2 - (P_x + Q_x) mod p, 
    y_3 = s(P_x - x_3) - P_y mod p */

    mpz_t s_square;
    mpz_init(s_square);
    mpz_mul(s_square, s, s);
    mpz_mod(s_square, s_square, ecc->prime);

    mpz_t temp1, temp2;
    mpz_inits(temp1, temp2, 0);

    /*x coordinates of the resulting point P + Q*/
    GF_addition(temp1, p->a, q->a, ecc->prime);
    GF_substruction(r->a, s_square, temp1, ecc->prime);

    /*y coordinates of the resulting point P + Q*/
    GF_substruction(temp2, p->a, r->a, ecc->prime);
    mpz_mul(temp2, temp2, s);
    mpz_mod(temp2, temp2, ecc->prime);
    GF_substruction(r->b, temp2, p->b, ecc->prime);

    mpz_clear(temp1);
    mpz_clear(temp2);
    mpz_clear(s);
    mpz_clear(s_square);
    mpz_clear(x_alpha);
    mpz_clear(y_alpha);
    mpz_clear(inverse_x_alpha);

    return(r);
};

/* This function if to perform point addition, where P = Q.

    s = (3(P_x)^2 + a) . (2P_y)^{-1}) mod p, if P = Q (point doubling)

    For furher information, 
    Christof Paar · Jan Pelzl
    Understanding Cryptography
    A Textbook for Students and Practitioners

 */

struct 
Point *point_doubling(struct Point *p, struct Elliptic_Curve *ecc){

     struct Point *r = create_point(); 

     /* Compute {P_x}^2 mod p */
     mpz_t px_square;
     mpz_init(px_square);

     mpz_mul(px_square, p->a, p->a);
     mpz_mod(px_square, px_square, ecc->prime);

     /* Compute 3 . {P_x}^2 mod p */
     mpz_t three;
     mpz_init_set_str(three, "3", 10);

     mpz_t px_square_3;
     mpz_init(px_square_3);
      
     mpz_mul(px_square_3, px_square, three);
     mpz_mod(px_square_3, px_square_3, ecc->prime);
     
     /* Compute 3 . {P_x}^2 + a mod p */
     mpz_t temp1;
     mpz_init(temp1);

     GF_addition(temp1, px_square_3, ecc->a, ecc->prime);

     /* Compute 2 . {P_y} mod p */
     mpz_t two;
     mpz_init_set_str(two, "2", 10); 

     mpz_t py_2;
     mpz_init(py_2);

     mpz_mul(py_2, p->b, two);
     mpz_mod(py_2, py_2, ecc->prime);

      /* Compute 2 . {P_x} mod p */ 
     mpz_t px_2;
     mpz_init(px_2);

     mpz_mul(px_2, p->a, two);
     mpz_mod(px_2, px_2, ecc->prime);

     /* Compute {2. {P_y}}^{-1} mod p */
     mpz_t inverse_py_2;
     mpz_init(inverse_py_2);
     mpz_invert(inverse_py_2, py_2, ecc->prime);
    
     /* Compute the slope, (3(P_x)^2 + a) . (2P_y)^{-1}) mod p */
     mpz_t s;
     mpz_init(s);

     mpz_mul(s, temp1, inverse_py_2);
     mpz_mod(s, s, ecc->prime);

     /* Compute s^2 mod q*/
     mpz_t s_square;
     mpz_init(s_square);
     mpz_mul(s_square, s, s);
     mpz_mod(s_square, s_square, ecc->prime);

     /*It is time to find x_3 and y_3.
    
     x_3 = s^2 - 2P_x mod p, 
     y_3 = s(P_x - x_3) - P_y mod p */

     GF_substruction(r->a, s_square, px_2, ecc->prime);

     mpz_t temp3;
     mpz_init(temp3);
     GF_substruction(temp3, p->a, r->a, ecc->prime);

     mpz_mul(temp3, s, temp3);
     mpz_mod(temp3, temp3, ecc->prime);

     GF_substruction(r->b, temp3, p->b, ecc->prime);

     return(r);

     mpz_clear(temp3);
     mpz_clear(px_square);
     mpz_clear(s_square);
     mpz_clear(s);
     mpz_clear(px_2);
     mpz_clear(py_2);
     mpz_clear(inverse_py_2);
     mpz_clear(temp1);
     mpz_clear(px_square_3);
     mpz_clear(three);
     mpz_clear(two);
    
};

/*This function is to perform scalar multiplication in ECC over prime finite field.
Given is an elliptic curve E, a generator point G = (G_x,G_y), an integer d such that 1 ≤ d ≤ #E 
with #E is the number of elements in elliptic curve E, we get a point T = (T_x,T_y) on the elliptic 
curve E such that dG = T.

                        G + G + G + ... + G = T
                        ----d times----
                this is nothing but a repeated addition 

p = the point we want to multiply.
dlp = how many times the point p will be multiplied
ecc = the elliptic curve that the point p is lying (prime will be used)
*/

struct 
Point *scalar_multiplication(struct Point *p, mpz_t dlp, struct Elliptic_Curve *ecc){

     /*Let us create a point, this will be return as a resulting point.*/
     struct Point *result = create_point();

     /*Let us copy the point p, which we want to multiply. */
     struct Point *p_copy = copy_point(p);

     struct Point *tmp1, *tmp2;

     /* Convert dlp to a string of digits in base 2. If str is NULL, the result string is allocated using the current allocation function. */
     char *val_str = mpz_get_str(NULL, 2, dlp);
     size_t len = strlen(val_str);

     int i;
     for (i = len - 1; i >= 0; i--) {
            
            /*If the processed bit is equal to 1, then we perform point addition.*/
		        if (val_str[i] == '1'){
            
			      tmp1 = point_addition(p_copy, result, ecc);
			      clear_point(result);
			      result = tmp1;
		        }

            /*Point doubling is default.*/
		        tmp2 = point_doubling(p_copy, ecc);
		        clear_point(p_copy);
		        p_copy = tmp2;
	}

	clear_point(p_copy);
	free(val_str);

	return result;

};

/*This function is to solve the modular equation y^2 = rhs (mod p) using the Shanks-Tonelli algorihm.*/

int 
square_root_mod_prime(mpz_t q, mpz_t rhs, struct Elliptic_Curve *ecc){

    /*Case when p|rhs, so rhs = 0 mod (p). Square root of 0 is 0.*/

    if(mpz_divisible_p(rhs, ecc->prime)) {         
      mpz_set_ui(q, 0);                          
      return 1;                                  
    }                                                            

    /*Case when rhs is not a quadratic residue. Return error, i.e, 0.                  */

    if(mpz_legendre(rhs, ecc->prime) != 1)         
        return 0;                                    

    /*Case when p = 3 (mod 4). If so, q = rhs ^ ((p+1) / 4) (mod p).                   */
  
    if(mpz_tstbit(ecc->prime, 1) == 1) {         
        mpz_set(q, ecc->prime);
        mpz_add_ui(q, q, 1);
        mpz_fdiv_q_2exp(q, q, 2);
        mpz_powm(q, rhs, q, ecc->prime);          
        return 1;
    }

    /* By factoring out powers of 2, find Q and S such that p - 1 = Q . 2^s with Q odd. */
    mpz_t z;
    mpz_init(z);
    unsigned int s;
   
    mpz_set(q, ecc->prime);                      
    mpz_sub_ui(q, q, 1);                        
    s = 0;                                       

    while(mpz_tstbit(q, s) == 0) s++;            
    mpz_fdiv_q_2exp(q, q, s);                      

    /*Search for a z in Z/pZ which is a quadratic non-residue. */
    mpz_set_ui(z, 2);                            
    while(mpz_legendre(z, ecc->prime) != -1)     
        mpz_add_ui(z, z, 1);                     
    
    /* z = z^q (mod p) */
    mpz_powm(z, z, q, ecc->prime);        

    /* q = (q+1)        
       q = (q+1) / 2   
       q = rhs^q (mod p) */

    mpz_add_ui(q, q, 1);                         
    mpz_fdiv_q_2exp(q, q, 1);                    
    mpz_powm(q, rhs, q, ecc->prime);              

    mpz_t inverse_of_rhs;
    mpz_init(inverse_of_rhs);

    mpz_t y;
    mpz_init(y);
    
    int i;

    for(;;) {
        mpz_powm_ui(y, q, 2, ecc->prime);        
        mpz_mul(y, y, inverse_of_rhs);
        mpz_mod(y, y, ecc->prime);               
        i = 0;
        /*As long as y is not equal to 1, */
        while(mpz_cmp_ui(y, 1) != 0) {
            i++;
            mpz_powm_ui(y, y, 2, ecc->prime);    
        }

        if(i == 0){                             
            return 1;
        }

        if(s-i == 1){                           
            mpz_mul(q, q, z);                  
        } 
        
        else{
           mpz_powm_ui(y, z, 1 << (s-i-1), ecc->prime);
            mpz_mul(q, q, y);
        }

        mpz_mod(q, q, ecc->prime);             
    }
}

void 
koblitz_embedding(mpz_t m, mpz_t k, struct Elliptic_Curve *ecc){

  unsigned long int k_new;
  k_new = mpz_get_ui(k);

  for(unsigned long int j = 0; j < k_new; j++){

    /* m . k */
    mpz_t m_multiplied_k;
    mpz_init(m_multiplied_k);
    mpz_mul(m_multiplied_k, m, k);

    /* x = m . k + j */
    mpz_t x;
    mpz_init(x);
    mpz_add_ui(x, m_multiplied_k, j);

    /* x^3 mod p */
    mpz_t cube_power;
    mpz_init(cube_power);
    mpz_powm_ui(cube_power,x, 3, ecc->prime);

    /* a . x */
    mpz_t a_multiplication_x;
    mpz_init(a_multiplication_x);
    mpz_mul(a_multiplication_x, x, ecc->a);

    /* a . x + b */
    mpz_t ax_addition_b;
    mpz_init(ax_addition_b);
    mpz_add(ax_addition_b, a_multiplication_x, ecc->b);

    /* x^3 + a . x + b */
    mpz_t rhs;
    mpz_init(rhs);
    GF_addition(rhs, cube_power, ax_addition_b, ecc->prime);
    
    int legendre;
    legendre = mpz_legendre(rhs, ecc->prime);

    if (legendre == 1){

      gmp_printf("\nFor x = m . k + j, x: %Zd\n", x);
      printf("y^2 = x^3 + x.a + b mod (p) is a quadratic residue.");

      mpz_t square;
      mpz_init(square);
      square_root_mod_prime(square, rhs, ecc);
      gmp_printf("\nSquare root of y^2 in mod (p) equals to: %Zd\n", square);

      gmp_printf("\nx coordinate of the message as an elliptic curve is: %Zd\n", x);
      gmp_printf("\ny coordinate of the message as an elliptic curve is: %Zd\n", square);

      /*I did not want to write all possible points. Whenever the program finds a succesfull point, it breaks.*/

      break;

      printf("----------------------------------------------------------------------------------------------------");

    } 
  }
}

void 
random_integer(mpz_t random_number, struct Elliptic_Curve *ecc){

  char *val_str = mpz_get_str(NULL, 2, ecc->n);
  size_t len = strlen(val_str);

  clock_t time = clock();
  gmp_randstate_t grt;
  gmp_randinit_default(grt);
  gmp_randseed_ui(grt, time);
  mpz_urandomb (random_number, grt, len);

  mpz_mod(random_number, random_number, ecc->n);

}

int
main() 
{
  
  /*Question 2.*/
  
  struct Elliptic_Curve *curve = (struct Elliptic_Curve*)malloc(sizeof(*curve));

  /*Elliptic Curve*/

  mpz_init_set_str(curve->a, "-3", 10); 
  mpz_init_set_str(curve->prime, "115792089210356248762697446949407573530086143415290314195533631308867097853951", 10);
  mpz_init_set_str(curve->b, "41058363725152142129326129780047268409114441015993725554835256314039467401291", 10);  
  mpz_init_set_str(curve->n, "115792089210356248762697446949407573529996955224135760342422259061068512044369", 10);

  gmp_printf("Elliptic curve component of prime: %Zd\n\n", curve->prime);
  gmp_printf("Elliptic curve component of a: %Zd\n\n", curve->a);
  gmp_printf("Elliptic curve component of b: %Zd\n\n", curve->b);
  gmp_printf("Elliptic curve component of n: %Zd\n\n", curve->n);

  printf("\n----------------------------------------------------------------------------------------------------");

  /*Generator G = (Gx, Gy)*/
  
  struct Point *gxgy = (struct Point*)malloc(sizeof(*gxgy));
  
	mpz_init_set_str(gxgy->a, "48439561293906451759052585252797914202762949526041747995844080717082404635286",10); 
	mpz_init_set_str(gxgy->b, "36134250956749795798585127919587881956611106672985015071877198253568414405109",10);

	gmp_printf("\nx coordinate of the pair G = (gx,gy): %Zd\n\n", gxgy->a);
  gmp_printf("\ny coordinate of the pair G = (gx,gy): %Zd\n\n", gxgy->b);

  printf("\n----------------------------------------------------------------------------------------------------");

  mpz_t randomAlice;
  mpz_init(randomAlice);
  random_integer(randomAlice, curve);
  gmp_printf("\nSecret key of Alice, k_prA: %Zd\n\n", randomAlice);

  struct Point *Alice_public_key = create_point(); 
  Alice_public_key = scalar_multiplication(gxgy, randomAlice, curve);

  gmp_printf("\nx coordinate of the Alice's public key: %Zd\n\n", Alice_public_key->a);
  gmp_printf("\ny coordinate of the Alice's public key: %Zd\n\n", Alice_public_key->b);

  printf("\n----------------------------------------------------------------------------------------------------");

  mpz_t randomBob;
  mpz_init(randomBob);
  random_integer(randomBob, curve);
  gmp_printf("\nSecret key of Bob, k_prB: %Zd\n\n", randomBob);

  struct Point *Bob_public_key = create_point(); 
  Bob_public_key = scalar_multiplication(gxgy, randomBob, curve);

  gmp_printf("\nx coordinate of the Bob's public key: %Zd\n\n", Bob_public_key->a);
  gmp_printf("\ny coordinate of the Bob's public key: %Zd\n\n", Bob_public_key->b);

  printf("\n----------------------------------------------------------------------------------------------------");

  /*Time to generate the common key.*/ 

  struct Point *what_Bob_generates = (struct Point*)malloc(sizeof(*what_Bob_generates));

  what_Bob_generates = scalar_multiplication(Alice_public_key, randomBob, curve);

  gmp_printf("\nx coordinate of what Bob generates: %Zd\n\n", what_Bob_generates->a);
  gmp_printf("\ny coordinate of what Bob generates: %Zd\n\n", what_Bob_generates->b);

  printf("\n----------------------------------------------------------------------------------------------------");

  struct Point *what_Alice_generates = (struct Point*)malloc(sizeof(*what_Alice_generates));

  what_Alice_generates = scalar_multiplication(Bob_public_key, randomAlice, curve);

  gmp_printf("\nx coordinate of what Alice generates: %Zd\n\n", what_Alice_generates->a);
  gmp_printf("\ny coordinate of what Alice generates: %Zd\n\n", what_Alice_generates->b);

  mpz_clear(randomAlice);
  mpz_clear(randomBob);
  clear_point(gxgy);
  clear_point(what_Bob_generates);
  clear_point(what_Alice_generates);
  clear_point(Bob_public_key);
  clear_point(Alice_public_key);

  printf("\n----------------------------------------------------------------------------------------------------");

  /*Question 1. */

  mpz_t m1, m2, m3;
  mpz_inits(m1, m2, m3, 0);
  
  EncodeFileIntoNumber("input_file", m1);
  gmp_printf("\nencoded message 1: %Zd\n", m1);

  printf("\n----------------------------------------------------------------------------------------------------");

  EncodeFileIntoNumber("input_file2", m2);
  gmp_printf("\nencoded message 2: %Zd\n", m2);

  printf("\n----------------------------------------------------------------------------------------------------");

  EncodeFileIntoNumber("input_file3", m3);
  gmp_printf("\nencoded number message 3: %Zd\n", m3);
  
  mpz_t k;
  mpz_init_set_str(k, "30",10);
  
  printf("\n----------------------------------------------------------------------------------------------------");

  koblitz_embedding(m1, k, curve);

  printf("\n----------------------------------------------------------------------------------------------------");

  koblitz_embedding(m2, k, curve);

  printf("\n----------------------------------------------------------------------------------------------------");

  koblitz_embedding(m3, k, curve);

  mpz_clear(m1);
  mpz_clear(m2);  
  mpz_clear(m3);  
  return 0;
}



