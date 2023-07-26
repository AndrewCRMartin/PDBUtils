/*************************************************************************

   Program:    protsurf
   File:       protsurf.c
   
   Version:    V1.2
   Date:       29.03.00
   Function:   Create contour plot of protein surface.
   
   Copyright:  (c) UCL, Dr. Andrew C. R. Martin 1994-2000
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   Phone:      (Home) +44 (0)1372 275775
   EMail:      andrew@stagleys.demon.co.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! 

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

   Compile with
   cc -o protsurf protsurf.c contour.c graphics.c -lbiop -lgen -lm

**************************************************************************

   Revision History:
   =================
   V1.0  17.06.94 Original
   V1.1  20.06.94 Added colour option; changed default grid & contours
   V1.2  29.03.00 Tidied up return value

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/macros.h"
#include "bioplib/pdb.h"

#include "contour.h"

/************************************************************************/
/* Defines and macros
*/
#define GRIDSTEP 2
#define MAXRAD   2.0
#define NCONT    5

/************************************************************************/
/* Globals
*/
int gColourPlot = 0;    /* Flag for producing colour plots              */

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *pdbfile, float *GridStep,
                  float *ContourStep);
void FindXYZLimits(PDB *pdb, 
                   float *xmin, float *xmax, float *ymin, float *ymax, 
                   float *zmin, float *zmax);
void FillGrid(float *Grid, int xsize, int ysize, float xmin, float ymin,
              float GridStep, PDB *pdb);
float CalcZ(float x, float y, PDB *p);
PDB *OpenAndReadPDB(char *file, FILE *Msgfp);
void Usage(void);
void PrintGrid(float *Grid, int xsize, int ysize);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for creating contour plots of protein surfaces.

   14.06.94 Original    By: ACRM
   17.06.94 Extended all limits by 4A
   29.03.00 Returns 0
*/
int main(int argc, char **argv)
{
   float GridStep    = GRIDSTEP,
         ContourStep = NCONT,
         xmin, xmax,
         ymin, ymax,
         zmin, zmax,
         *Grid;
   int   xsize, 
         ysize;
   char  pdbfile[160];
   PDB   *pdb,
         *p;
   
   if(ParseCmdLine(argc, argv, pdbfile, &GridStep, &ContourStep))
   {
      if((pdb = OpenAndReadPDB(pdbfile, stderr)) != NULL)
      {
         FindXYZLimits(pdb,&xmin,&xmax,&ymin,&ymax,&zmin,&zmax);

         /* Move the protein along z, so all z-coordinates are +ve      */
         for(p=pdb; p!=NULL; NEXT(p))
            p->z -= zmin;
         zmax -= zmin;
         zmin  = 0;

         /* Extend all limits by 4A each way                            */
         xmax += 4.0;
         ymax += 4.0;
         xmin -= 4.0;
         ymin -= 4.0;
         
         xsize = 1 + (int)((xmax-xmin)/GridStep);
         ysize = 1 + (int)((ymax-ymin)/GridStep);

         if((Grid = (float *)malloc(xsize * ysize * sizeof(float)))!=NULL)
         {
            FillGrid(Grid,xsize,ysize,xmin,ymin,GridStep,pdb);
            psOpen();
            Contour(Grid,xsize,ysize,-ContourStep);
#ifdef DEBUG
            PrintGrid(Grid,xsize,ysize);
#endif
            psClose();
         }
         else
         {
            fprintf(stderr,"No memory for grid\n");
         }
      }
      else
      {
         fprintf(stderr,"No atoms read from PDB file\n");
      }
   }
   else
   {
      Usage();
   }

   return(0);
}

/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *pdbfile, 
                     float *GridStep, float *ContourStep)
   -------------------------------------------------------
   Input:   int   argc          Argument count
            char  **argv        Argument array
   Output:  char  *pdbfile      PDB file for input
            float *GridStep     Grid step size
            float *ContourStep  +ve Number of contours
                                -ve Contour step size
   Returns: BOOL                Command line OK?

   Parses the command line

   14.06.94 Original    By: ACRM
   20.06.94 Sets gColourPlot
*/
BOOL ParseCmdLine(int argc, char **argv, char *pdbfile, float *GridStep,
                  float *ContourStep)
{
   argc--;
   argv++;
   
   if(argc < 1)
      return(FALSE);
   
   while(argc > 1)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'g':
            argc--;
            argv++;
            if((sscanf(argv[0],"%f",GridStep))==0)
               return(FALSE);
            break;
         case 'c':
            argc--;
            argv++;
            if((sscanf(argv[0],"%f",ContourStep))==0)
               return(FALSE);
            break;
         case 'm':
            gColourPlot = 1;
            break;
         default:
            return(FALSE);
         }
      }
      else
      {
         return(FALSE);
      }
      
      argc--;
      argv++;
   }
   
   strcpy(pdbfile,argv[0]);
   
   return(TRUE);
}

/************************************************************************/
/*>void FindXYZLimits(PDB *pdb, 
                      float *xmin, float *xmax, 
                      float *ymin, float *ymax, 
                      float *zmin, float *zmax)
   --------------------------------------------
   Input:   PDB    *pdb       PDB linked list
   Output:  float  *xmin      Limits of box including centres of atoms
                   *xmax
                   *ymin
                   *ymax
                   *zmin
                   *zmax

   Find the min and max coordinates of a PDB linked list. The limits
   returned are those of the atoms centres, no account is taken of
   atom radii.

   14.06.94 Original    By: ACRM
*/
void FindXYZLimits(PDB *pdb, 
                   float *xmin, float *xmax, float *ymin, float *ymax, 
                   float *zmin, float *zmax)
{
   PDB *p;

   *xmin = *xmax = pdb->x;
   *ymin = *ymax = pdb->y;
   *zmin = *zmax = pdb->z;

   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(p->x < *xmin) *xmin = p->x;
      if(p->y < *ymin) *ymin = p->y;
      if(p->z < *zmin) *zmin = p->z;

      if(p->x > *xmax) *xmax = p->x;
      if(p->y > *ymax) *ymax = p->y;
      if(p->z > *zmax) *zmax = p->z;
   }
}
         
/************************************************************************/
/*>void FillGrid(float *Grid, int xsize, int ysize, float xmin, 
                 float ymin, float GridStep, PDB *pdb)
   ------------------------------------------------------------
   Input:   int    xsize     X dimension of grid
            int    ysize     Y dimension of grid
            float  xmin      Min x value
            float  ymin      Min y value
            float  GridStep  Step size to take across the grid
            PDB    *pdb      PDB linked list
   Output:  float  *Grid     The z-value grid arranged in a 1D array

   Fill in the grid with z values calculated from the PDB linked list.

   14.06.94 Original    By: ACRM
*/
void FillGrid(float *Grid, int xsize, int ysize, float xmin, float ymin,
              float GridStep, PDB *pdb)
{
   PDB *p;
   float x, y, z;
   int  i, j;

   for(j=0; j<ysize; j++)
   {
      y = ymin + j*GridStep;
      
      for(i=0; i<xsize; i++)
      {
         x = xmin + i*GridStep;
         
         Grid[(ysize-j-1)*xsize + i] = (float)(0.0);


         for(p=pdb; p!=NULL; NEXT(p))
         {
            if(p->x >= x-MAXRAD && p->x <= x+MAXRAD &&
               p->y >= y-MAXRAD && p->y <= y+MAXRAD)
            {
               z = CalcZ(x,y,p);
               
               if(z > Grid[(ysize-j-1)*xsize + i])
                  Grid[(ysize-j-1)*xsize + i] = z;
            }
         }
      }
   }
}

/************************************************************************/
/*>float CalcZ(float x, float y, PDB *p)
   -------------------------------------
   Input:   float   x      X-grid coordinate
            float   y      Y-grid coordinate
            PDB     *p     PDB structure pointer
   Returns: float          Z value at grid point.

   Calculate a z value given a grid point and PDB pointer. If the grid
   point is outside the atom's boundaries, returns the value 0.0

   14.06.94 Original    By: ACRM
*/
float CalcZ(float x, float y, PDB *p)
{
   float rad = 1.7,
        xoff, yoff, zoff,
        z   = 0.0;

   /* Correct the radius for oxygen                                     */
   if(p->atnam[0] == 'O') rad = 1.35;
   
   /* Check the point is within the circle                              */
   if((((x - p->x)*(x - p->x)) + ((y - p->y)*(y - p->y))) <= rad*rad)
   {
      xoff = x - p->x;
      yoff = y - p->y;
      
      zoff = (float)sqrt((double)(rad*rad - xoff*xoff - yoff*yoff));
      z = xoff + p->z;
   }
   
   return(z);
}
   
/************************************************************************/
/*>PDB *OpenAndReadPDB(char *file, FILE *Msgfp)
   --------------------------------------------
   Opens a PDB file specified by name and reads it with ReadPDB. Error
   messages are sent to the specified file. If NULL, no messages will
   be issued. Returns a pointer to the PDB linked list or NULL if failed.

   10.10.93 Original    By: ACRM
   11.10.93 Added legal NULL Msgfp
*/
PDB *OpenAndReadPDB(char *file, FILE *Msgfp)
{
   FILE *fp = NULL;
   PDB  *pdb = NULL;
   int  natoms;

   if((fp=fopen(file,"r"))==NULL)
   {
      if(Msgfp!=NULL) 
         fprintf(Msgfp,"Unable to open file: %s\n",file);
   }
   else
   {
      pdb = ReadPDB(fp, &natoms);
      if(pdb == NULL || natoms == 0)
      {
         if(Msgfp!=NULL) 
            fprintf(Msgfp,"No atoms read from file: %s\n",file);
         if(pdb!=NULL) FREELIST(pdb,PDB);
         pdb = NULL;
      }
      fclose(fp);
   }
   return(pdb);
}

/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a usage message

   14.06.94 Original    By: ACRM
   17.06.94 Updated
   20.06.94 Added -m option
*/
void Usage(void)
{
   fprintf(stderr,"\nProtSurf V1.1 (c) 1994 Andrew C.R. Martin, UCL\n");
   fprintf(stderr,"Uses contour code from DDJ. ProtSurf is freely \
distributable providing\n");
   fprintf(stderr,"no profit is made; the contour code is usable under \
the conditions of DDJ\n\n");
   fprintf(stderr,"Usage: protsurf [-g <gridstep>] [-c <ncontour>] [-m] \
<file.pdb>\n");
   fprintf(stderr,"       -g Grid stepsize (Default: %4.1f)\n", 
           (double)GRIDSTEP);
   fprintf(stderr,"       -c Number of contours (Default: %4.1f)\n", 
           (double)NCONT);
   fprintf(stderr,"          Use a negative number of contours to \
specify contour separation\n");
   fprintf(stderr,"       -m Multi-colour plot\n\n");
   fprintf(stderr,"Generates PostScript output of a protein surface \
contour plot.\n\n");
}

/************************************************************************/
/*>void PrintGrid(float *Grid, int xsize, int ysize)
   -------------------------------------------------
   Input:   float   *Grid    Array of grid points
            int     xsize    X-dimension of grid
            int     ysize    Y-dimension of grid

   Prints the grid for debugging

   14.06.94 Original    By: ACRM
*/
void PrintGrid(float *Grid, int xsize, int ysize)
{
   int x, y;
   
   for(x=0; x<xsize; x++)
   {
      for(y=0; y<ysize; y++)
      {
         fprintf(stderr,"%8.3f ",Grid[y*xsize+x]);
      }
      fprintf(stderr,"\n");
   }
}

