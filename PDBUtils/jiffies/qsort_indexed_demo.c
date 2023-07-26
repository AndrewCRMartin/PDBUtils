#include <stdlib.h>
#include <stdio.h>

int D[] = {5,2,9,3,1,7};

int compare(const void *a, const void *b)
{
   if(D[*(int *)a] > D[*(int *)b])
      return(1);
   else if(D[*(int *)a] == D[*(int *)b])
      return(0);
   else
      return(-1);
}


int main(int argc, char **argv)
{
   int I[] = {0,1,2,3,4,5};
   int i;
   
   qsort((void *)I, (size_t)6, (size_t)sizeof(int), &compare);
   
   for(i=0; i<6; i++)
   {
      printf("%d\n", D[I[i]]);
   }
   return(0);
}

