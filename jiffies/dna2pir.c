#include <stdio.h>
#include <string.h>
#include "bioplib/seq.h"
#include "bioplib/macros.h"

int main(int argc, char **argv)
{
   FILE *in;
   char dna[4096],
        buffer[256];
   int  i=0;
   
   dna[0] = '\0';
   
   if((in = fopen(argv[1],"r"))!=NULL)
   {
      printf(">P1;%s\n",argv[1]);
      printf("Translated by - dna2pir\n");
      
      while(fgets(buffer,256,in))
      {
         TERMINATE(buffer);
         strcat(dna, buffer);
      }
   }

   for(i=0; i<strlen(dna); i+=3)
   {
      printf("%c",DNAtoAA(dna+i));
      if(i && !(i%90))
         printf("\n");
   }
   if(!(i%90))
      printf("\n");
   
   return(0);
}

      
      
   
   
