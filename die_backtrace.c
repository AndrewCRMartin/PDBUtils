#include <stdio.h>
#include <stdlib.h>

void die(void);
int main(int argc, char **argv);
int backtrace(void **, int);
char **backtrace_symbols(void **, int);


int main(int argc, char **argv)
{
   die();
   return(0);
}

void die(void)
{
   void *stack[20];
   char **functions;
   int count, i;
   
   count=backtrace(stack, 20);
   functions = backtrace_symbols(stack, count);
   for(i=0; i<count; i++)
   {
      fprintf(stderr,"Frame %2d: %s\n", i, functions[i]);
   }
   free(functions);
   exit(1);
}
