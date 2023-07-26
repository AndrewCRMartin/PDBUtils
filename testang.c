/*************************************************************************

   Program:    testang
   File:       
   
   Version:    
   Date:       
   Function:   
   
   Copyright:  (c) Dr. Andrew C. R. Martin 1995
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   Phone:      (Home) +44 (0)1372 275775
               (Work) +44 (0)171 387 7050 X 3284
   EMail:      INTERNET: martin@biochem.ucl.ac.uk
               
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
   Program to find possible range of CA-CA-CA angles by spinning
   around omega1, phi2, psi2, omega2

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/seq.h"
#include "bioplib/angle.h"
#include "bioplib/macros.h"

/************************************************************************/
/* Defines
*/
#define MAXTOR 20   /* Max number of torsions (and therefore residues)  */

/* Parameters from CHARMm                                               */
#define ANGLE_C   117.5*PI/180.0
#define ANGLE_N   120.0*PI/180.0
#define ANGLE_CA  111.6*PI/180.0
#define DIST_CA_C 1.52
#define DIST_C_N  1.33
#define DIST_N_CA 1.49

/************************************************************************/
/* Macros
*/

/************************************************************************/
/* Type definitions
*/
struct _entry
{
   struct _entry *next,
                 *prev;
   REAL          x,
                 y,
                 z;
   char          atnam[8];
}  ;
typedef struct _entry ENTRY;

/************************************************************************/
/* Globals
*/
FILE *gHitfp = NULL,
     *gTorfp = NULL,
     *gOutfp = NULL;

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ProcessCmdLine(int argc, char **argv, BOOL *calphas, char *hitfile,
                    char *torfile, char *outfile, BOOL *DoOmega,
                    BOOL *SmallOut, BOOL *DoNNFile);
BOOL OpenFiles(char *hitfile, char *torfile, char *outfile);
void DoProcessing(BOOL calphas, BOOL DoOmega, BOOL SmallOut, 
                  BOOL DoNNFile);
int GetTorsions(FILE *fp, int nres, REAL *phi, REAL *psi, BOOL DoOmega,
                REAL *omega);
ENTRY *BuildFragment(int ntor, REAL *phi, REAL *psi, REAL *omega);
BOOL BuildAtom(ENTRY *p, ENTRY *q, ENTRY *r, REAL theta, REAL bond, 
               REAL phi, REAL *coords);
void WriteFragment(FILE *gOutfp, ENTRY *entry, BOOL calphas, 
                   BOOL SmallOut, char *seq);
char *GetSequence(char *buffer);
void StoreAngle(BOOL Display, ENTRY *entry);



/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for building loops from torsion data produced by searchdb

   08.11.93 Original   By: ACRM
   10.11.93 Modified to support omega angle
   07.06.94 Modified to support neural nets files
*/
int main(int argc, char **argv)
{
   REAL  omega1,
         omega2,
         phi,
         psi,
         step = (REAL)5.0*PI/180.0,
         Phi[8],
         Psi[8],
         Omega[8];
   ENTRY *entry;
   int   i;

   for(i=0;i<8;i++)
   {
      Phi[i]   = PI;
      Psi[i]   = PI;
      Omega[i] = PI;
   }

   for(omega1 = 0; omega1 <= PI; omega1 += PI)
   {
/*      Omega[0] = omega1;  */
      
      for(omega2 = 0; omega2 <= PI; omega2 += PI)
      {
         Omega[1] = omega2; 
         
         for(phi = -PI; phi <= PI; phi += step)
         {
            Phi[1] = phi;

            for(psi = -PI; psi <= PI; psi += step)
            {
               Psi[1] = psi;
               if((entry = BuildFragment(3, Phi, Psi, Omega))==NULL)
               {
                  fprintf(stderr,"No memory for entry\n");
                  return(1);
               }
               else
               {
                  StoreAngle(FALSE,entry);
                  FREELIST(entry,ENTRY);
               }
            }
         }
         
      }
   }
   
   StoreAngle(TRUE,NULL);

   return(0);
}


/************************************************************************/
/*>ENTRY *BuildFragment(int ntor, REAL *phi, REAL *psi, REAL *omega)
   -----------------------------------------------------------------
   Builds a fragment from a set of torsion angles. Assuming a fixed
   trans omega angle for the peptide bond:
   each phi angle defines the coordinates of this residues's C (and O)
   each psi angle defines the coordinates of the next residue's N
   each omega angle defines the coordinates of the next residue's CA
   The data are built into a ENTRY linked list.

   09.11.93 Original   By: ACRM
   10.11.93 Added omega parameter
*/
ENTRY *BuildFragment(int ntor, REAL *phi, REAL *psi, REAL *omega)
{
   ENTRY *entry, *p;
   REAL  coords[3];
   int   i;

   /* Check input                                                       */
   if(ntor == 0) return(NULL);

   /* Initialise an ENTRY structure for the start N                     */
   INITPREV(entry,ENTRY);
   if(entry == NULL) return(NULL);

   /* Fill in the data for the start N                                  */
   p = entry;
   strcpy(p->atnam,"N   ");
   p->x = (REAL)(-DIST_C_N);
   p->y = (REAL)0.0;
   p->z = (REAL)0.0;

   /* And for the C-alpha                                               */
   ALLOCNEXTPREV(p,ENTRY);
   if(p == NULL)
   {
      FREELIST(entry,ENTRY);
      return(NULL);
   }
   strcpy(p->atnam,"CA  ");
   p->x = (REAL)0.0;
   p->y = (REAL)0.0;
   p->z = (REAL)0.0;
   
   /* And for the C                                                     */
   ALLOCNEXTPREV(p,ENTRY);
   if(p == NULL)
   {
      FREELIST(entry,ENTRY);
      return(NULL);
   }
   strcpy(p->atnam,"C   ");
   p->x = (REAL)(DIST_CA_C * sin((double)(PI*ANGLE_CA/180.0)));
   p->y = (REAL)(DIST_CA_C * cos((double)(PI*ANGLE_CA/180.0)));
   p->z = (REAL)0.0;
   
   for(i=0; i<ntor; i++)
   {
      if(i!=0)
      {
         D("Built phi C atom\n");
         /* Build the C position using the phi data                     */
         if(BuildAtom(p->prev->prev, p->prev, p,   /* C, N, CA coords   */
                      ANGLE_CA, DIST_CA_C,         /* Constants         */
                      phi[i],                      /* Torsion angle     */
                      coords))                     /* Results           */
	 {
            ALLOCNEXTPREV(p,ENTRY);
            if(p == NULL)
            {
               FREELIST(entry,ENTRY);
               return(NULL);
            }
            strcpy(p->atnam,"C   ");
            p->x = coords[0];
            p->y = coords[1];
            p->z = coords[2];
	 }
      }

      if(i != (ntor-1))
      {
         /* Build the N position using the psi data                     */
         if(BuildAtom(p->prev->prev, p->prev, p,   /* N, CA, C coords   */
                      ANGLE_C, DIST_C_N,           /* Constants         */
                      psi[i],                      /* Torsion angle     */
                      coords))                     /* Results           */
	 {
            D("Built psi N atom\n");
            ALLOCNEXTPREV(p,ENTRY);
            if(p == NULL)
            {
               FREELIST(entry,ENTRY);
               return(NULL);
            }
            strcpy(p->atnam,"N   ");
            p->x = coords[0];
            p->y = coords[1];
            p->z = coords[2];
	 }

         /* Build the CA position using a 180 degree omega angle        */
         if(BuildAtom(p->prev->prev, p->prev, p,   /* CA, C, N coords   */
                      ANGLE_N, DIST_N_CA,          /* Constants         */
                      omega[i],                    /* Torsion angle     */
                      coords))                     /* Results           */
	 {
            D("Built omega CA atom\n");
            ALLOCNEXTPREV(p,ENTRY);
            if(p == NULL)
            {
               FREELIST(entry,ENTRY);
               return(NULL);
            }
            strcpy(p->atnam,"CA  ");
            p->x = coords[0];
            p->y = coords[1];
            p->z = coords[2];
	 }
      }
   }
   return(entry);
}

/************************************************************************/
/*>BOOL BuildAtom(ENTRY *p, ENTRY *q, ENTRY *r, REAL ang, REAL dist, 
                  REAL tor, REAL coords[3])
   -----------------------------------------------------------------
   Build coords for an atom given the coords of 3 anticedants, an angle,
   distance and torsion. Based on code from CARTX2 in CONGEN

   09.11.93 Original   By: ACRM
*/
#define ETA 1.0e-7
#define ETA2 ETA*ETA


BOOL BuildAtom(ENTRY *p, ENTRY *q, ENTRY *r, REAL theta, REAL bond, 
               REAL phi, REAL *coords)
{
   REAL stht, ctht,
        sphi, cphi,
        bsin,
        x1, y1, z1,
        x2, y2, z2,
        x3, y3, z3,
        x4, y4, z4,
        lyz1, ovlyz1,
        yy4, zz4,
        y1o, z1o,
        lxz22, l2, lxz2,
        ovl2, ovlxz2,
        x2o, z2o, xz2o, y2o, xx1, xx4;

   if(p==NULL || q==NULL || r==NULL) return(FALSE);

   stht  = (REAL)sin((double)(PI-theta));
   ctht  = (REAL)cos((double)(PI-theta));
   sphi  = (REAL)sin((double)phi);
   cphi  = (REAL)cos((double)phi);
   bsin  = bond * stht;

   x4    = bond * ctht;
   y4    = bsin * cphi;
   z4    = bsin * sphi;

   x3    = r->x;
   y3    = r->y;
   z3    = r->z;

   x1    = p->x - x3;
   y1    = p->y - y3;
   z1    = p->z - z3;

   x2    = q->x - x3;
   y2    = q->y - y3;
   z2    = q->z - z3;

   lxz22 = x2*x2 + z2*z2;
   l2    = (REAL)sqrt((double)(lxz22+y2*y2));
   lxz2  = (REAL)sqrt((double)lxz22);

   if(l2 < ETA)
   {
      fprintf(stderr,"Atoms 2 & 3 too close!\n");
      ovl2 = (REAL)1.0/ETA;
   }
   else
   {
      ovl2 = (REAL)1.0/l2;
   }

   if(lxz2 < ETA)
   {
      xx1 = x1;
      x2o = (REAL)1.0;
      z2o = (REAL)0.0;
   }
   else
   {
      ovlxz2 = (REAL)1.0/lxz2;
      x2o    = x2 * ovlxz2;
      z2o    = z2 * ovlxz2;
      xx1    = x1*x2o + z1*z2o;
      z1     = z1*x2o - x1*z2o;
   }

   xz2o   = lxz2 * ovl2;
   y2o    = y2   * ovl2;

   x1     = -xx1*xz2o - y1*y2o;
   y1     = xx1*y2o   - y1*xz2o;

   lyz1   = (REAL)sqrt((double)(y1*y1 + z1*z1));
   ovlyz1 = (REAL)1.0 / lyz1;

   y1o    = y1 * ovlyz1;
   z1o    = z1 * ovlyz1;

   yy4    = y1o*y4 - z1o*z4;
   zz4    = y1o*z4 + z1o*y4;
   xx4    = y2o*yy4 - xz2o*x4;

   y4     = -xz2o*yy4 - y2o*x4;
   x4     = x2o*xx4 - z2o*zz4;
   z4     = z2o*xx4 + x2o*zz4;

   coords[0] = x4 + x3;
   coords[1] = y4 + y3;
   coords[2] = z4 + z3;

   return(TRUE);
}

/************************************************************************/
void StoreAngle(BOOL Display, ENTRY *entry)
{
   ENTRY *CAs[3], *e;
   int   CACount = 0;
   REAL  ang;
   static REAL MinAng = 1000.0,
               MaxAng = -1000.0;
   
   if(Display)
   {
      printf("Min angle is %f (%f)\n",MinAng,180.0*MinAng/PI);
      printf("Max angle is %f (%f)\n",MaxAng,180.0*MaxAng/PI);
   }
   else
   {
      for(e=entry; e!=NULL; NEXT(e))
      {
         if(!strncmp(e->atnam,"CA  ",4))
         {
            CAs[CACount++] = e;
            if(CACount==3)break;
         }
      }
      
      ang = angle(CAs[0]->x, CAs[0]->y, CAs[0]->z,
                  CAs[1]->x, CAs[1]->y, CAs[1]->z,
                  CAs[2]->x, CAs[2]->y, CAs[2]->z);
      if(ang > MaxAng)
         MaxAng = ang;
      if(ang < MinAng)
         MinAng = ang;
   }
}

