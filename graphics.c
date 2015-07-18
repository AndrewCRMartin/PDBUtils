/* PostScript graphics for Contour from DDJ June 92 p 95
*/

#include <stdio.h>
#include "bioplib/SysDefs.h"

typedef struct
{
   float x,y;
}  LIST;

void Polyline(int inn, LIST *inlist)
{
   LIST *list = inlist;
   int  n     = inn;
   
   if(n<2)  return;
   
   /* First draw a filled area                                          */
   printf("newpath\n");
   printf("%.6f %.6f moveto\n",list->x, 1.0-list->y);
   list++;

   while(--n)
   {
      printf("%.6f %.6f lineto\n",list->x, 1.0-list->y);
      list++;
   }
   printf("closepath\n");
   printf("fill\n");
   printf("stroke\n\n");

   /* Now redraw the path over the top                                  */
   list = inlist;
   n    = inn;
   
   printf("0.0 setgray\n");
   printf("newpath\n");
   printf("%.6f %.6f moveto\n",list->x, 1.0-list->y);
   list++;

   while(--n)
   {
      printf("%.6f %.6f lineto\n",list->x, 1.0-list->y);
      list++;
   }
   printf("closepath\n");
   printf("stroke\n\n");
}


void ContourText(char *s, float x, float y)
{
   printf("0.0 setgray\n");
   printf("%.6f %.6f moveto (%s) show\n",x,1.0-y,s);
}

void psOpen(void)
{
   printf("%%!\n");
   printf("save\n\n");
   printf("/Helvetica findfont 0.015 scalefont setfont\n\n");

   printf("72 252 translate\n");
   printf("468 468 scale\n");
   printf("0.001 setlinewidth\n\n");
   printf("newpath\n");
   printf("0 0 moveto\n");
   printf("0 1 lineto\n");
   printf("1 1 lineto\n");
   printf("1 0 lineto\n");
   printf("closepath\n");
   printf("stroke\n");
   printf("clippath\n\n");
/*   printf("0.00001 setlinewidth\n\n");
*/
   printf("0.001 setlinewidth\n\n");
}

void psClose(void)
{
   printf("restore\n");
   printf("showpage\n");
}

