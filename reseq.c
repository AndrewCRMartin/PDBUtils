/*************************************************************************

   Program:    reseq
   File:       reseq.c
   
   Version:    V1.1
   Date:       08.07.96
   Function:   Change the sequence of a PDB file by MOP
   
   Copyright:  (c) Dr. Andrew C. R. Martin 1996
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   Phone:      (Home) +44 (0)1372 275775
               (Work) +44 (0)171 419 3890
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
   V1.0  06.02.96 Original   By: ACRM
   V1.1  08.07.96 Prints error message from RepSChain()

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/pdb.h"
#include "bioplib/seq.h"
#include "bioplib/general.h"

/************************************************************************/
/* Defines and macros
*/
#define CHITAB   "chitab.dat"
#define REFCOORD "coor"

#define MAXSEQS 16

#define MAXBUFF 160


/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL DoReseq(FILE *seq_fp, FILE *pdb_in, FILE *pdb_out);
char *BuildSeqString(char **seqs, int nchain);
BOOL ParseCmdLine(int argc, char **argv, char *seqfile, char *infile, 
                  char *outfile);
void Usage(void);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for resequencing

   06.02.96 Original   By: ACRM
*/
int main(int argc, char **argv)
{
   FILE *seq_fp,
        *pdb_in  = stdin,
        *pdb_out = stdout;
   char SeqFile[MAXBUFF],
        InFile[MAXBUFF],
        OutFile[MAXBUFF];

   if(ParseCmdLine(argc, argv, SeqFile, InFile, OutFile))
   {
      if(OpenStdFiles(InFile, OutFile, &pdb_in, &pdb_out))
      {
         if((seq_fp=fopen(SeqFile,"r"))==NULL)
         {
            fprintf(stderr,"Can't read sequence file: %s\n",SeqFile);
            return(1);
         }

         if(!DoReseq(seq_fp, pdb_in, pdb_out))
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
/*>BOOL DoReseq(FILE *seq_fp, FILE *pdb_in, FILE *pdb_out)
   -------------------------------------------------------
   Reads the files and calls the resequencing code

   06.02.96 Original   By: ACRM
   08.07.96 Prints the error message from RepSChain()
*/
BOOL DoReseq(FILE *seq_fp, FILE *pdb_in, FILE *pdb_out)
{
   char *seqs[MAXSEQS],
        *seq;
   BOOL punct, error;
   int  nchain, natom;
   PDB  *pdb;

   if((nchain=ReadPIR(seq_fp,FALSE,seqs,MAXSEQS,NULL,&punct,&error))==0)
   {
      fprintf(stderr,"Can't read sequence from PIR file\n");
      return(FALSE);
   }

   if((pdb=ReadPDB(pdb_in, &natom))==NULL)
   {
      fprintf(stderr,"Can't read atoms from PDB file\n");
      return(FALSE);
   }

   if((seq = BuildSeqString(seqs, nchain))==NULL)
   {
      fprintf(stderr,"No memory for sequence string\n");
      return(FALSE);
   }
   
   if(!RepSChain(pdb, seq, CHITAB, REFCOORD))
   {
      fprintf(stderr,"Resequencing failed: %s\n",gRSCError);
      return(FALSE);
   }

   WritePDB(pdb_out, pdb);
   
   return(TRUE);
}


/************************************************************************/
/*>char *BuildSeqString(char **seqs, int nchain)
   ---------------------------------------------
   Takes a multi-chain sequence specification and builds it into a single
   string. Returns a malloc'd character pointer

   05.02.96 Original    By: ACRM
*/
char *BuildSeqString(char **seqs, int nchain)
{
   int  i, 
        seqlen;
   char *string;
   
   /* Find the total length of the string                               */
   seqlen = 0;
   for(i=0; i<nchain; i++)
      seqlen += strlen(seqs[i]);

   /* Allocate this much space                                          */
   if((string = (char *)malloc(seqlen * sizeof(char)))==NULL)
      return(NULL);
   
   /* Copy the strings into the buffer                                  */
   string[0] = '\0';
   for(i=0;i<nchain;i++)
      strcat(string,seqs[i]);

   /* Return the buffer pointer                                         */
   return(string);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *seqfile, char *infile, 
                     char *outfile)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *seqfile     PIR sequence file
            char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
   Returns: BOOL                Success?

   Parse the command line
   
   05.02.96 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *seqfile, char *infile, 
                  char *outfile)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   
   if(!argc)
      return(FALSE);

   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are 1, 2 or 3 arguments left               */
         if(argc > 3 || argc < 1)
            return(FALSE);
         
         /* Copy the first to seqfile                                   */
         strcpy(seqfile, argv[0]);


         /* If there's another, copy it to infile                       */
         argc--;
         argv++;
         if(argc)
            strcpy(infile, argv[0]);

         /* If there's another, copy it to outfile                      */
         if(argc)
         {
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

   06.02.96 Original   By: ACRM
*/
void Usage(void)
{      
   fprintf(stderr,"\nreseq V1.1 (c) 1996, Dr. Andrew C.R. Martin, UCL\n");

   fprintf(stderr,"\nUsage: reseq <seq.pir> [<in.pdb> [<out.pdb>]]\n");

   fprintf(stderr,"\nApplies the sequence from the PIR file to the input \
PDB file, writing\n");
   fprintf(stderr,"it to the output file. If filenames are not specified \
for the PDB files,\n");
   fprintf(stderr,"standard input/output are used.\n\n");
}

