/*************************************************************************

   Program:    pdbreseq
   File:       pdbreseq.c
   
   Version:    V1.0
   Date:       11.11.94
   Function:   Resequence a PDB file
   
   Copyright:  (c) Dr. Andrew C. R. Martin 1994
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   Phone:      (Home) +44 (0372) 275775
   EMail:      INTERNET: amartin@scitec.adsp.sub.org
                         martin@bsm.bioc.ucl.ac.uk
               UUCP:     ...{uunet|rutgers}!cbmehq!cbmuk!scitec!amartin
               JANET:    martin@uk.ac.ucl.bioc.bsm
               
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
#include <string.h>
#include <stdlib.h>
#include "bioplib/pdb.h"
#include "bioplib/general.h"

/************************************************************************/
/* Defines and macros
*/
/*
#define CHITAB    "/home/bsm/martin/data/chilink"
#define REFCOORD  "/home/bsm/martin/data/coor"
*/
#define CHITAB    "chilink"
#define REFCOORD  "coor"

#define MAXBUFF 160

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
BOOL ParseCmdLine(int argc, char **argv, char **sequence, char *chitab, 
                  char *refcoord, char *infile, char *outfile);
void Usage(void);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for protein resquencing filter

   11.11.94 Original   By: ACRM
*/
int main(int argc, char **argv)
{
   char infile[MAXBUFF],
        outfile[MAXBUFF],
        chitab[MAXBUFF],
        refcoord[MAXBUFF],
        *sequence = NULL;
   FILE *in  = stdin,
        *out = stdout;
   PDB  *pdb;
   int  natoms;

   if(ParseCmdLine(argc, argv, &sequence, chitab, refcoord, 
                   infile, outfile))
   {
      if(OpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdb = ReadPDB(in,&natoms)) != NULL)
         {
            if(RepSChain(pdb,sequence,chitab,refcoord))
            {
               WritePDB(out, pdb);
            }
            else
            {
               fprintf(stderr,"Resequencing failed: %s\n",gRSCError);
               return(1);
            }
         }
         else
         {
            fprintf(stderr,"No atoms read from PDB file.\n");
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
/*>BOOL ParseCmdLine(int argc, char **argv, char **sequence, 
                     char *chitab, char *refcoord, char *infile, 
                     char *outfile)
   -------------------------------------------------------------
   Parse the command line

   11.11.94 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char **sequence, char *chitab, 
                  char *refcoord, char *infile, char *outfile)
{
   argc--;
   argv++;

   if(!argc)
      return(FALSE);

   infile[0]  = outfile[0] = '\0';
   strcpy(chitab,CHITAB);
   strcpy(refcoord,REFCOORD);
   *sequence = NULL;
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'c':
            argc--;
            argv++;
            strncpy(chitab,argv[0],MAXBUFF);
            chitab[MAXBUFF-1] = '\0';
            break;
         case 'r':
            argc--;
            argv++;
            strncpy(refcoord,argv[0],MAXBUFF);
            refcoord[MAXBUFF-1] = '\0';
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that a sequence has been specified                    */
         if(argc < 1)
            return(FALSE);

         /* Allocate space and copy the sequence                        */
         *sequence = (char *)malloc((strlen(argv[0])+1)*sizeof(char));
         if(*sequence == NULL)
         {
            fprintf(stderr,"Unable to allocate memory for sequence\n");
            exit(1);
         }
         strcpy(*sequence,argv[0]);

         argc--;
         argv++;

         if(argc)
         {
            /* Check that there are only 1 or 2 arguments left          */
            if(argc > 2)
               return(FALSE);
            
            /* Copy the first to infile                                 */
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
   
   return(TRUE);
}

/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a usage message

   11.11.94 Original    By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\nPDBReseq V1.0 (c) 1994, Dr. Andrew C.R. Martin, \
UCL.\n");
   fprintf(stderr,"Freely distributable if no profit is made in so \
doing.\n\n");
   fprintf(stderr,"Usage: pdbreseq [-c chitab] [-r refcoords] sequence \
[in.pdb] [out.pdb]\n");
   fprintf(stderr,"       -c Specify equivalent chi table\n");
   fprintf(stderr,"       -r Specify reference coordinate set\n");
   fprintf(stderr,"I/O is through stdin/stdout if files are not \
specified.\n\n");
}
