/*************************************************************************

   Program:    pdbfit
   File:       pdbfit.c
   
   Version:    V1.0
   Date:       12.12.01
   Function:   Very crude simple fit of 2 identical structures
   
   Copyright:  (c) University of Reading / Dr. Andrew C. R. Martin 2001
   Author:     Dr. Andrew C. R. Martin
   Address:    School of Animal and Microbial Sciences,
               The University of Reading,
               Whiteknights,
               P.O. Box 228,
               Reading RG6 6AJ.
               England.
   Phone:      +44 (0)118 987 5123 Extn. 7022
   Fax:        +44 (0)118 931 0180
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
   PDB  *pdb1, *pdb2;
   REAL rm[3][3];
   
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
   
   if(natoms1 != natoms2)
   {
      fprintf(stderr,"Non-identical PDB lists\n");
      return(1);
   }
   
   if(!FitNCaCOPDB(pdb1, pdb2, rm))
   {
      fprintf(stderr,"Fitting failed\n");
      return(1);
   }
   
   WritePDB(stdout,pdb2);
   
   return(0);
}

void Usage(void)
{
   fprintf(stderr,"\npdbfit V1.0 12.12.01 ACRM\n");
   fprintf(stderr,"Usage: pdbfit file1.pdb file2.pdb [>file_out.pdb]\n");
   fprintf(stderr,"Simple crude program to fit two identical (content) \
PDB files\n\n");
}
