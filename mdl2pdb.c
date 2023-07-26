/*************************************************************************

   Program:    mdl2pdb
   File:       mdl2pdb.c
   
   Version:    V1.1
   Date:       26.07.95
   Function:   Convert Charmm MDL to PDB
   
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
   Simple C program to convert CHARMM >MDL format to PDB.
   N.B. Does not handle chain names!

**************************************************************************

   Usage:
   ======
   Compile with:

   cc -o mdl2pdb mdl2pdb.c

**************************************************************************

   Revision History:
   =================
   V1.0  21.07.95 Original    By: ACRM
   V1.1  26.07.95 Added return to main()

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
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
void DoProcessing(FILE *in, FILE *out);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile);
void Usage(void);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program to MDL to PDB conversion

   21.07.95 Original    By: ACRM
*/
int main(int argc, char **argv)
{
   FILE   *in   = stdin,
          *out  = stdout;
   char   infile[MAXBUFF],
          outfile[MAXBUFF];

   if(ParseCmdLine(argc, argv, infile, outfile))
   {
      if(OpenStdFiles(infile, outfile, &in, &out))
      {
         DoProcessing(in, out);
         return(0);
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
/*>void DoProcessing(FILE *in, FILE *out)
   --------------------------------------
   Does the actual work of MDL to PDB conversion

   21.07.95 Original    By: ACRM
*/
void DoProcessing(FILE *in, FILE *out)
{   
   BOOL First = TRUE;
   int  atnum,
        resnum,
        ijunk;
   REAL x, y, z, fjunk;
   char resnam[16],
        atnam[16],
        segid[16];
   char buffer[MAXBUFF];

   while(fgets(buffer, MAXBUFF, in))
   {
      if(strchr(buffer,'*'))
      {
         fprintf(out,"REMARK   1 %s",buffer);
      }
      else if(First)
      {
         First = FALSE;
      }
      else
      {
         sscanf(buffer,"%d %d %s %s %lf %lf %lf %s %d %lf",
                &atnum,
                &resnum,
                resnam,
                atnam,
                &x,
                &y,
                &z,
                segid,
                &ijunk,
                &fjunk);
         fprintf(out,"ATOM  %5d  %-4s%-4s %4d    %8.3f%8.3f%8.3f  1.00 20.00\n",
                 atnum,
                 atnam,
                 resnam,
                 resnum,
                 x,
                 y,
                 z);
      }
   }
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
   
   21.07.95 Original    By: ACRM
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
/*>void Usage(void)
   ----------------
   Prints a usage message

   21.07.95 Original    By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\nmdl2pdb V1.0 (c) Dr. Andrew C.R. Martin, UCL. \
Freely distributable\n");

   fprintf(stderr,"\nUsage: mdl2pdb [file.mdl [file.pdb]]\n");

   fprintf(stderr,"\nConverts a Charmm .mdl (or MDL .mol) file to PDB \
format.\n");
   fprintf(stderr,"N.B. Ignores chain identifiers\n\n");
}

