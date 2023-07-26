/*************************************************************************

   Program:    calcpdbdist
   File:       calcpdbdist.c
   
   Version:    V1.0
   Date:       11.12.17
   Function:   Calculate the distance between C-alpha atoms of the 
               two residues specified on the command line
   
   Copyright:  (c) UCL / Dr. Andrew C. R. Martin 2017
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Institute of Structural & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   EMail:      andrew.martin@ucl.ac.uk
               andrew@bioinf.org.uk
               
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

   Should this be modified to work only with sidechain atoms???

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0   11.12.18   Original   By: ACRM

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <math.h>
#include "bioplib/macros.h"
#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/pdb.h"
#include "bioplib/general.h"
#include "bioplib/array.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *res1, char *res2,
                  char *infile, char *outfile);
void Usage(void);
REAL CalcDistance(PDB *pdb, char *res1, char *res2);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**
   Main program

-  11.12.17   Original   By: ACRM
*/
int main(int argc, char **argv)
{
   FILE *in  = stdin,
        *out = stdout;
   char InFile[MAXBUFF],
        OutFile[MAXBUFF],
        res1[MAXBUFF],
        res2[MAXBUFF];

   if(ParseCmdLine(argc, argv, res1, res2, InFile, OutFile))
   {
      if(blOpenStdFiles(InFile, OutFile, &in, &out))
      {
         PDB *pdb;
         int natom;
         
         if((pdb=blReadPDB(in, &natom))!=NULL)
         {
            REAL d;
            
            d = CalcDistance(pdb, res1, res2);
            if(d < 0.0)
            {
               return(1);
            }
            fprintf(out, "%.3f\n", d);
         }
         else
         {
            fprintf(stderr,"calcpdbdistances: (Error) No atoms read from PDB \
file\n");
            return(1);
         }
      }
      else
      {
         fprintf(stderr,"calcpdbdistances: (Error) No atoms read from PDB \
file\n");
            return(1);
         
      }
      
   }
   else
   {
      Usage();
   }
   
   return(0);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *res1, char *res2,
                     char *infile, char *outfile)
   -----------------------------------------------------------------
*//**

   \param[in]      argc                 Argument count
   \param[in]      **argv               Argument array
   \param[out]     *res1                One key residue
   \param[out]     *res2                Second key residue
   \param[out]     *infile              Input file (or blank string)
   \param[out]     *outfile             Output file (or blank string)
   \return                              Success?

   Parse the command line

-  11.12.17  Original   By: ACRM   
*/
BOOL ParseCmdLine(int argc, char **argv, char *res1, char *res2, 
                  char *infile, char *outfile)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   res1[0]                = '\0';
   res2[0]                = '\0';
   
   if(!argc)
   {
      return(FALSE);
   }
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         if (argv[0][2]!='\0')
         {
           return(FALSE);
         }
         else
         {            
            switch(argv[0][1])
            {
            case 'h':
               return(FALSE);
               break;
            default:
               return(FALSE);
               break;
            }
         }         
      }
      else
      {
         /* Check that there are 2, 3 or 4 arguments left               */
         if(argc < 2 || argc > 4)
            return(FALSE);
         
         /* Copy the first to res1 and second one to res2               */
         strcpy(res1, argv[0]);
         argc--;
         argv++;
         strcpy(res2, argv[0]);
         argc--;
         argv++;
         
         /* Copy the first to infile                                    */
         if(argc)
         {
            strcpy(infile, argv[0]);
            argc--;
            argv++;
         }
         
         /* If there's another, copy it to outfile                      */
         if(argc)
         {
            strcpy(outfile, argv[0]);
            argc--;
            argv++;
         }
         
         return(TRUE);
      }
      
      argc--;
      argv++;
   }
   return(TRUE);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
*//**
   Prints a usage message

-  11.12.17  Original   By: ACRM
*/
void Usage(void)
{
   fprintf(stderr, "\ncalcpdbdistance V1.0 (c) 2017 UCL, Dr. Andrew C.R. \
Martin\n");

   fprintf(stderr, "\nUsage: calcpdbdistance res1 res2 \
[in.pdb [out.txt]]\n");

   fprintf(stderr, "Calculates the distance between the C-alpha atoms \
of the two specified\n");
   fprintf(stderr, "residues. I/O is through standard input / output \
if files are not\n");
   fprintf(stderr, "specified.\n");

   fprintf(stderr, "\n");
   blPrintResSpecHelp(stderr);
   fprintf(stderr, "\n");
}

/************************************************************************/
/*>REAL CalcDistance(PDB *pdb, char *res1, char *res2)
   ---------------------------------------------------
*//**
   \param[in]    *out    Output file pointer
   \param[in]    *pdb    PDB linked list
   \param[in]    res1    First residue spec
   \param[in]    res2    Second residue spec
   \return               Distance (<0 on failure)

   Calculates the distance between the C-alphas of the two 
   specified residues.

-  11.12.17   Original   By: ACRM
*/
REAL CalcDistance(PDB *pdb, char *res1, char *res2)
{
   PDB *r1 = NULL,
       *r2 = NULL;
   REAL d  = 0.0;
   
   
   if((r1 = blFindResidueSpec(pdb, res1))==NULL)
   {
      fprintf(stderr,"Error (calcpdbdistance): residue not found (%s)\n",
              res1);
      return(-1.0);
   }
   if((r2 = blFindResidueSpec(pdb, res2))==NULL)
   {
      fprintf(stderr,"Error (calcpdbdistance): residue not found (%s)\n",
              res2);
      return(-1.0);
   }

   if((r1 = blFindAtomInRes(r1, "CA"))==NULL)
   {
      fprintf(stderr,"Error (calcpdbdistance): CA atom not found in %s\n",
              res1);
      return(-1.0);
   }
   if((r2 = blFindAtomInRes(r2, "CA"))==NULL)
   {
      fprintf(stderr,"Error (calcpdbdistance): CA atom not found in %s\n",
              res2);
      return(-1.0);
   }

   d = DIST(r1, r2);
   return(d);
}

