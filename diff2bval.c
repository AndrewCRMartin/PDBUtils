/*************************************************************************

   Program:    diff2bval
   File:       diff2bval.c
   
   Version:    V1.0
   Date:       21.08.95
   Function:   Write atomic distances between 2 PDB files into the BVal
               column. Writes both input files to the output file.
   
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

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================

*************************************************************************/
/* Includes
*/

/************************************************************************/
/* Defines and macros
*/
#include <stdio.h>
#include <math.h>
#include "bioplib/pdb.h"
#include "bioplib/macros.h"

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   (Badly written) whole program for writing atom distances into B-val
   column.

   21.08.95 Original    By: ACRM
*/
int main(int argc, char **argv)
{
   PDB  *pdb1, 
        *pdb2,
        *p, *q;
   FILE *in1, 
        *in2, 
        *out;
   int  natoms;
   
   if(argc!=4)
   {
      fprintf(stderr,"\ndiff2bval (c) 1995, Dr. Andrew C.R. Martin, \
UCL\n");
      fprintf(stderr,"Usage: diff2bval <in1.pdb> <in2.pdb> <out.pdb>\n");

      fprintf(stderr,"\nTakes 2 PDB files of same structure which are \
superimposed and calculates\n");
      fprintf(stderr,"distance between equivalent atoms writing the \
difference in the BVal\n");
      fprintf(stderr,"column of the output file which contains both \
structures.\n\n");

      return(1);
   }

   if((in1=fopen(argv[1],"r"))==NULL)
   {
      fprintf(stderr,"Unable to open input file %s\n",argv[1]);
      return(1);
   }
   if((in2=fopen(argv[2],"r"))==NULL)
   {
      fprintf(stderr,"Unable to open input file %s\n",argv[2]);
      return(1);
   }
   if((out=fopen(argv[3],"w"))==NULL)
   {
      fprintf(stderr,"Unable to open output file %s\n",argv[3]);
      return(1);
   }

   pdb1 = ReadPDB(in1, &natoms);
   pdb2 = ReadPDB(in2, &natoms);

   for(p=pdb1, q=pdb2; p!=NULL && q!=NULL; NEXT(p), NEXT(q))
      p->bval = q->bval = DIST(p,q);

   WritePDB(out,pdb1);
   WritePDB(out,pdb2);

   return(0);
}

   
   
   
   
