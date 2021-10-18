#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>

///*  linux  */ #include <gmp.h>
/* windows */ #include <mpir.h>

/* 
*  Change 1 to 0 in
*  #define PRINT_TEXT_BLOCK 1
*  below if you don't want to print inside the procedures. */
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

int
main() 
{
  mpz_t A, B;
  mpz_inits(A, B, 0);
  
  EncodeFileIntoNumber("input_file", A);
  gmp_printf("\nencoded number A: %Zd\n", A);
  
  // Do something
  mpz_add_ui(B, A, 1);
  gmp_printf("B: %Zd\n", B);
  mpz_sub_ui(B, B, 1);
  gmp_printf("B: %Zd\n\n", B);
  
  DecodeNumberIntoFile("output_file", B);
  
  //mpz_clears(A, B, 0);
}
