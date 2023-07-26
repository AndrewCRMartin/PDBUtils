/*************************************************************************

   Program:    getresidue
   File:       getresidue.c
   
   Version:    V1.0
   Date:       14.07.11
   Function:   Extract the residue name for a specifed ID
   
   Copyright:  (c) UCL / Dr. Andrew C. R. Martin 2011
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
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
   V1.0  14.07.11 Original   By: ACRM

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include "bioplib/pdb.h"
#include "bioplib/general.h"
#include "bioplib/macros.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF  160

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
void Usage(void);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  char *resspec);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for extracting selected chains from a PDB file

   07.02.97 Original   By: ACRM
   06.04.09 Added lowercase option
   22.05.09 Added keepHeader
   29.06.09 Added atomsOnly
*/
int main(int argc, char **argv)
{
   char InFile[MAXBUFF],
        OutFile[MAXBUFF],
        resspec[MAXBUFF];
   FILE *in  = stdin,
        *out = stdout;
   PDB  *pdb = NULL;
   int  natoms;
   
   if(ParseCmdLine(argc, argv, InFile, OutFile, resspec))
   {
      if(blOpenStdFiles(InFile, OutFile, &in, &out))
      {
         if((pdb=blReadPDB(in, &natoms)) == NULL)
         {
            fprintf(stderr,"No atoms read from input PDB file\n");
            return(1);
         }
         else
         {
            PDB *p;
            p = blFindResidueSpec(pdb, resspec);
            if(p!=NULL)
            {
               fprintf(out, "%s\n", p->resnam);
            }
            else
            {
               fprintf(out, "NULL\n");
            }
         }
      }
   }
   else
   {
      Usage();
   }

   return(0);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a usage message

   14.07.11 Original   By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\ngetresidue V1.0 (c) 2011 Dr. Andrew C.R. Martin, UCL\n");

   fprintf(stderr,"\nUsage: getresidue [chain]resnum[insert] \
[in.pdb [outfile]]\n");

   fprintf(stderr,"\nGetresidue identifies the amino acid at a specified \
position in a\n");
   fprintf(stderr,"PDB file.\n\n");
} 


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                     char *resspec)
   ----------------------------------------------------------------------
   Input:   int    argc        Argument count
            char   **argv      Argument array
   Output:  char   *infile     Input filename (or blank string)
            char   *outfile    Output filename (or blank string)
            char   *resspec    Residue ID
   Returns: BOOL               Success

   Parse the command line

   14.07.11 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  char *resspec)
{
   argc--;
   argv++;
   
   infile[0] = outfile[0] = resspec[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
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
      else
      {
         /* Check that there are 1-3 arguments left                     */
         if(argc > 3)
            return(FALSE);
         
         /* Copy the first to resspec                                   */
         strcpy(resspec, argv[0]);
         
         /* If there's another, copy it to infile                       */
         argc--;
         argv++;
         if(argc)
         {
            strcpy(infile, argv[0]);
         
            /* If there's another, copy it to outfile                   */
            argc--;
            argv++;
            if(argc)
               strcpy(outfile, argv[0]);
         }

         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   if(!resspec[0])
      return(FALSE);

   return(TRUE);
}

