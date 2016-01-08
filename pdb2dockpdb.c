/*************************************************************************

   Program:    pdb2dockpdb
   File:       pdb2dockpdb.c
   
   Version:    V1.1
   Date:       31.05.02
   Function:   Convert standard PDB to Dock extended format
   
   Copyright:  (c) Dr. Andrew C. R. Martin 1996-2002
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
   V1.0  14.02.96 Original
   V1.1  31.05.02 Changed PDB field from 'junk' to 'record_type'

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/general.h"

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
void Usage(void);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile);
void WriteDockPDB(FILE *fp,
                  PDB  *pdb);
void WriteDockPDBRecord(FILE *fp,
                        PDB  *pdb);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   14.02.96 Original   By: ACRM
*/
int main(int argc, char **argv)
{
   FILE *in  = stdin,
        *out = stdout;
   char infile[MAXBUFF],
        outfile[MAXBUFF];
   PDB  *pdb;
   int  natoms;

   if(ParseCmdLine(argc,argv,infile,outfile))
   {
      if(OpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdb = ReadPDB(in, &natoms))!=NULL)
            WriteDockPDB(out, pdb);
         else
            return(1);
      }
      else
      {
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
/*>void Usage(void)
   ----------------
   14.02.96 Original   By: ACRM
   31.05.02 V1.1
*/
void Usage(void)
{
   fprintf(stderr,"\npdb2dockpdb V1.1 (c) 1996-2002, \
Dr. Andrew C.R. Martin, UCL\n");

   fprintf(stderr,"\nUsage: pdb2dockpdb [in.pdb [out.dpdb]]\n");

   fprintf(stderr,"\nConverts a standard PDB file to Dock extended PDB \
format. Assumes that\n");
   fprintf(stderr,"charges have been placed in the BVal column. All \
atom types are set\n");
   fprintf(stderr,"to zero\n\n");
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
   Returns: BOOL                Success?

   Parse the command line
   
   14.02.96 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   
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
         /* Check that there are only 1 or 2 arguments left             */
         if(argc > 2)
            return(FALSE);
         
         /* Copy the first to infile                                    */
         strcpy(infile, argv[0]);
         
         /* If there's another, copy it to outfile                      */
         argc--;
         argv++;
         if(argc)
            strcpy(outfile, argv[0]);
            
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}

/************************************************************************/
/*>void WriteDockPDB(FILE *fp, PDB *pdb)
   -------------------------------------
   Input:   FILE *fp   PDB file pointer to be written
            PDB  *pdb  PDB linked list to write

   Write a PDB linked list in Dock exteneded PDB format by calls to 
   WriteDockPDBRecord()

   14.02.96 Original
*/
void WriteDockPDB(FILE *fp,
                  PDB  *pdb)
{
   PDB   *p;
   char  PrevChain[8];
   
   strcpy(PrevChain,pdb->chain);

   for(p = pdb ; p ; NEXT(p))
   {
      if(strncmp(PrevChain,p->chain,1))
      {
         /* Chain change, insert TER card                               */
         fprintf(fp,"TER   \n");
         strcpy(PrevChain,p->chain);
      }
      WriteDockPDBRecord(fp,p);
   }
   fprintf(fp,"TER   \n");
}

/************************************************************************/
/*>void WriteDockPDBRecord(FILE *fp, PDB *pdb)
   -------------------------------------------
   Input:   FILE  *fp     PDB file pointer to be written
            PDB   *pdb    PDB linked list record to write

   Write a Dock extended PDB record

   14.02.96 Original
*/
void WriteDockPDBRecord(FILE *fp,
                        PDB  *pdb)
{
   fprintf(fp,"%-6s%5d  %-4s%-4s%1s%4d%1s   \
%8.3f%8.3f%8.3f%8.3f%8.3f%3d\n",
           pdb->record_type,
           pdb->atnum,
           pdb->atnam,
           pdb->resnam,
           pdb->chain,
           pdb->resnum,
           pdb->insert,
           pdb->x,
           pdb->y,
           pdb->z,
           pdb->bval,
           pdb->occ,
           0);
}
