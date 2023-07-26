/*************************************************************************

   Program:    pdbcafit
   File:       pdbcafit.c
   
   Version:    V1.0
   Date:       12.05.10
   Function:   Very crude simple fit of 2 identical structures
   
   Copyright:  (c) UCL / Dr. Andrew C. R. Martin 2010
   Author:     Dr. Andrew C. R. Martin
   EMail:      andrew@bioinf.org.uk
               
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
#include <stdio.h>
#include <stdlib.h>
#include "bioplib/pdb.h"
#include "bioplib/macros.h"

/************************************************************************/
/* Defines and macros
*/

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
void Usage(void);

/************************************************************************/
int main(int argc, char **argv)
{
   FILE *in1, *in2;
   int  natoms1, natoms2;
   PDB  *pdb1, *pdb2, *p;
   REAL rm[3][3], rms;
   
   if(argc != 3)
   {
      Usage();
      return(1);
   }
   
   /* Open files */
   if((in1 = fopen(argv[1], "r"))==NULL)
   {
      fprintf(stderr,"Can't read %s\n", argv[1]);
      return(1);
   }
   if((in2 = fopen(argv[2], "r"))==NULL)
   {
      fprintf(stderr,"Can't read %s\n", argv[2]);
      return(1);
   }

   /* Read coordinates */
   if((pdb1=ReadPDB(in1, &natoms1))==NULL)
   {
      fprintf(stderr,"Can't read atoms from %s\n",argv[1]);
      return(1);
   }
   if((pdb2=ReadPDB(in2, &natoms2))==NULL)
   {
      fprintf(stderr,"Can't read atoms from %s\n",argv[2]);
      return(1);
   }

   pdb1 = SelectCaPDB(pdb1);
   pdb2 = SelectCaPDB(pdb2);
   
   natoms1 = natoms2 = 0;
   for(p=pdb1; p!=NULL; NEXT(p)) natoms1++;
   for(p=pdb2; p!=NULL; NEXT(p)) natoms2++;
   
   if(natoms1 != natoms2)
   {
      fprintf(stderr,"Non-identical PDB lists\n");
      return(1);
   }
   
   if(!FitCaPDB(pdb1, pdb2, rm))
   {
      fprintf(stderr,"Fitting failed\n");
      return(1);
   }
   
   rms = CalcRMSPDB(pdb1, pdb2);
   
   printf("%.3f\n", rms);
   

/*   WritePDB(stdout,pdb2); */
   
   return(0);
}

void Usage(void)
{
   fprintf(stderr,"\npdbcafit V1.0 12.05.10 ACRM\n");
   fprintf(stderr,"Usage: pdbcafit file1.pdb file2.pdb\n");
   fprintf(stderr,"Simple crude program to fit CAs of two identical \
(content) PDB files\n");
   fprintf(stderr,"Returns prints the CA-RMSD\n\n");
}
