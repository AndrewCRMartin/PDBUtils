#include <stdio.h>
#include <ctype.h>

int main(int argc, char **argv)
{
   int ch;
   
   if(argc!=1)
   {
      fprintf(stderr,"Usage: printable <in >out\n");
      return(1);
   }
   
   while(((ch=fgetc(stdin))!=EOF) && !feof(stdin))
   {
      if(isprint(ch))
      {
         fputc(ch, stdout);
      }
      else if(ch==(char)13)
      {
         fputc('\n', stdout);
      }
   }
   return(0);
}
