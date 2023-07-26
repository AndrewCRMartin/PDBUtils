/*************************************************************************

   Program:    pdbclash
   File:       pdbclash.c
   
   Version:    V1.2
   Date:       03.10.17
   Function:   Display a list of atom clashes from a PDB file
   
   Copyright:  (c) UCL / Dr. Andrew C. R. Martin 2019
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

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   12.01.19 V1.0    Original   By: ACRM

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/access.h"

/************************************************************************/
/* Defines and macros
*/
#define RESBUFF 16
#define MAXBUFF 1024
#define DEF_BINWIDTH 1.0
#define DEF_RADFILE "radii.dat"
#define DATA_ENV    "DATADIR"
#define ISBACKBONE(p) (!strncmp((p)->atnam, "N   ", 4) || \
                       !strncmp((p)->atnam, "CA  ", 4) || \
                       !strncmp((p)->atnam, "C   ", 4) || \
                       !strncmp((p)->atnam, "O   ", 4))

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  char *radFile, BOOL *checkBB, REAL *tol);
void Usage(void);
int main(int argc, char **argv);
BOOL Analyze(FILE *outfile, PDB *pdb, BOOL checkBB, REAL tol);
void CheckClashes(FILE *out,
                  PDB *res1, PDB *nextRes1,
                  PDB *res2, PDB *nextRes2,
                  BOOL checkBB, REAL tol);


/************************************************************************/
int main(int argc, char **argv)
{
   char infile[MAXBUFF],
      outfile[MAXBUFF],
      radFile[MAXBUFF];
   int  natoms;
   FILE *in  = stdin,
        *out = stdout,
        *fpRad = NULL;
   PDB  *pdb, *pdbin;
   BOOL noenv,
        checkBB = FALSE;
   REAL tol = (REAL)0.0;

   strncpy(radFile, DEF_RADFILE, MAXBUFF);
   
   if(ParseCmdLine(argc, argv, infile, outfile, radFile, &checkBB, &tol))
   {
      /* Open the radius file                                           */
      if((fpRad=blOpenFile(radFile, DATA_ENV, "r", &noenv))==NULL)
      {
         fprintf(stderr, "Error (pdbclash): Unable to open radius file, \
%s\n", radFile);
         if(noenv)
         {
            fprintf(stderr, "              Environment variable %s \
not set\n", DATA_ENV);
         }
         return(1);
      }

      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdbin = blReadPDB(in, &natoms)) == NULL)
         {
            fprintf(stderr, "Unable to read PDB file\n");
            return(1);
         }


         /* Remove waters                                               */
         pdb = blStripWatersPDBAsCopy(pdbin, &natoms);
         FREELIST(pdbin, PDB);
         
         /* Set the atom radii in the linked list                       */
         blSetAtomRadii(pdb, fpRad);
         
         /* Perform the analysis                                        */
         if(!Analyze(out, pdb, checkBB, tol))
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
BOOL Analyze(FILE *out, PDB *pdb, BOOL checkBB, REAL tol)
{
   PDB *res1, *res2, *nextRes1, *nextRes2;

#ifdef DEBUG
   PDB *p;
   for(p=pdb; p!=NULL; NEXT(p))
      p->bval = p->radius;
   blWritePDB(out, pdb);
#endif

   /* Step through one residue at a time                                */
   for(res1=pdb; res1!=NULL; res1=nextRes1)
   {
      /* Find the next residue                                          */
      nextRes1 = blFindNextResidue(res1);

      /* Step through subsequent residues one at a time                 */
      for(res2=nextRes1; res2!=NULL; res2=nextRes2)
      {
         /* Find the next residue                                       */
         nextRes2 = blFindNextResidue(res2);

         /* And check for clashes...                                    */
         CheckClashes(out, res1, nextRes1, res2, nextRes2, checkBB, tol);
      }
   }
   
   return(TRUE);
}

void CheckClashes(FILE *out,
                  PDB *res1, PDB *nextRes1,
                  PDB *res2, PDB *nextRes2,
                  BOOL checkBB, REAL tol)
{
   PDB *p, *q;
   
   for(p=res1; p!=nextRes1; NEXT(p))
   {
      for(q=res2; q!=nextRes2; NEXT(q))
      {
         REAL sumVDWR   = (p->radius + q->radius);
         REAL sumVDWRSq = (sumVDWR-tol) * (sumVDWR-tol);
         REAL distSq    = DISTSQ(p, q);

         /* If we are looking at sequence-adjacent residues             */
         if((res2 == nextRes1) &&
            CHAINMATCH(res2->chain, nextRes1->chain))
         {
            if(checkBB)
            {
               /* If we are looking at two backbone atoms then skip it  */
               if(ISBACKBONE(p) && ISBACKBONE(q))
                  continue;
            }
            else
            {
               continue;
            }
         }

         /* Ignore Cys-(SG/CB) - Cys-(SG/CB) as this is probably a disulphide if the
            distance is low
         */
         if(!strncmp(p->resnam, "CYS", 3) &&
            !strncmp(q->resnam, "CYS", 3))
         {
            if(!strncmp(p->atnam, "SG  ", 4) && !strncmp(q->atnam, "SG  ", 4))
               continue;
            if(!strncmp(p->atnam, "SG  ", 4) && !strncmp(q->atnam, "CB  ", 4))
               continue;
            if(!strncmp(p->atnam, "CB  ", 4) && !strncmp(q->atnam, "SG  ", 4))
               continue;
         }
         
         
         if(distSq < sumVDWRSq)
         {
            REAL dist = sqrt(distSq);
            char resspec1[16],
                 resspec2[16];
            blBuildResSpec(p, resspec1);
            blBuildResSpec(q, resspec2);
            
            fprintf(out, "%5s:%s clash with %5s:%s \
Dist: %.2f Allowed: %.2f Clash: %.2f\n",
                    resspec1, p->atnam,
                    resspec2, q->atnam,
                    dist, sumVDWR, (sumVDWR-dist));
         }
      }
   }
}

/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                     char *radFile, BOOL *checkBB, REAL *tol) 
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
            char   *radFile     Radius file
            BOOL   *checkBB     If we are looking at adjecent residues
                                should we also check that atoms are
                                both backbone
            REAL   *tol         Tolerance
   Returns: BOOL                Success?

   Parse the command line
   
   12.01.19 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  char *radFile, BOOL *checkBB, REAL *tol)
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
         case 'r':
            argc--;
            argv++;
            if(!argc)
               return(FALSE);
            strncpy(radFile, argv[0], MAXBUFF);
            break;
         case 't':
            argc--;
            argv++;
            if(!argc)
               return(FALSE);
            sscanf(argv[0], "%lf", tol);
            break;
         case 'b':
            *checkBB = TRUE;
            break;
         case 'h':
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are <2 arguments left                      */
         if(argc > 2)
            return(FALSE);
         
         if(argc)
         {
            /* Copy the next to infile                                  */
            strcpy(infile, argv[0]);
            argc--;
            argv++;
            
            /* If there's another, copy it to outfile                   */
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
*/
void Usage(void)
{
   printf("\npdbclash V1.0 (c) 2019 UCL, Andrew C.R. Martin\n");
   printf("\nUsage: pdbclash [-b][-r radii.dat][-t x] [input.pdb \
[output.txt]]\n");
   printf("         -b When skipping 'bad' contacts between adjacent \
residues,\n");
   printf("            should we only skip bad contacts between \
backbone atoms?\n");
   printf("         -t Specify tolerence for allowed clashes [0.0]\n");

   printf("\nTakes a PDB file and looks for clashes between atoms \
taking into account\n");
   printf("the VDW radii of the atoms.\n\n");
}

