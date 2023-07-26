#include <stdio.h>
#include <math.h>
/* #include <exec/types.h> */
#include "bioplib/SysDefs.h"

void Contour(float *data,
             int dim_x,
             int dim_y,
             double inc);


float data[400];

main(int argc, char **argv)
{
   float *s;
   int i,j;
   double x,y;
   
   s=data;
   for(i=0, s=data; i<20; i++)
   {
      for(j=0; j<20; j++)
      {
         x=((double)i - 9.5)/4.0;
         y=((double)j - 9.5)/4.0;
         
         *s++ = (float)(10.0*cos(x) * cos(y));
      }
   }
   psOpen();
   Contour(data,20,20,2.0);
   psClose();
}
