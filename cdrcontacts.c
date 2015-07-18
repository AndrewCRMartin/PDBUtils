/*************************************************************************

   Program:    cdrcontacts
   File:       cdrcontacts.c
   
   Version:    V1.1
   Date:       03.12.14
   Function:   Find contacts between CDRs and framework
   
   Copyright:  (c) UCL / Dr. Andrew C. R. Martin 2014
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
   V1.0   16.11.14  Original   By: ACRM
   V1.1   03.12.14  Output now includes residue names

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bioplib/macros.h"
#include "bioplib/pdb.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXLABEL      16
#define MAXBUFF       400
#define DEF_DISTCUT   (REAL)4.0
#define MAXBONDDISTSQ (REAL)5.29  /* i.e. 2.3A (Nx...O(x-1))            */

typedef struct _zone
{
   char cdr[MAXLABEL],
        start[MAXLABEL],
        stop[MAXLABEL];
} ZONE;

   

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
void RunAnalysis(FILE *out, PDB *pdb, ZONE *cdrs, REAL dist);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  REAL *dist, BOOL *doResList);
void Usage(void);
BOOL InCDR(PDB *p, ZONE *cdrs, char *theCDR);
REAL MakesContact(PDB *res1, PDB *res2, REAL distSq);
BOOL InPDBZoneSpec(PDB *p, char *start, char *stop);
void PrintResList(FILE *out, PDB *pdb);


/************************************************************************/
int main(int argc, char **argv)
{
   /* Initialize the zones for CDRs and the whole of an antibody        */
   ZONE cdrs[] = {{"CDR-L1", "L24" , "L34" },
                  {"CDR-L2", "L50" , "L56" },
                  {"CDR-L3", "L89" , "L97" },
                  {"CDR-H1", "H26" , "H35D"},
                  {"CDR-H2", "H50" , "H65" },
                  {"CDR-H3", "H95" , "H102"},
                  {""      , ""    , ""    }};
   FILE *in  = stdin,
        *out = stdout;
   char infile[MAXBUFF],
        outfile[MAXBUFF];
   REAL dist = DEF_DISTCUT;
   BOOL doResList = FALSE;
   

   if(ParseCmdLine(argc, argv, infile, outfile, &dist, &doResList))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         PDB *pdb;
         int natoms;
         
         if((pdb=blReadPDBAtoms(in, &natoms))!=NULL)
         {
            RunAnalysis(out, pdb, cdrs, dist);
            if(doResList)
               PrintResList(out, pdb);
            FREELIST(pdb, PDB);
         }
      }
      else
      {
         fprintf(stderr,"Error: Unable to open I/O files\n");
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
/*>void PrintResList(FILE *out, PDB *pdb)
   --------------------------------------
   04.12.14 Original   By: ACRM
*/
void PrintResList(FILE *out, PDB *pdb)
{
   PDB *p;
   
   fprintf(out, "\n\nRESLIST\n");
   for(p=pdb; p!=NULL; p=blFindNextResidue(p))
   {
      fprintf(out, "%s%d%s %s\n",
              p->chain, p->resnum, p->insert, p->resnam);
   }
}


/************************************************************************/
/*>void RunAnalysis(FILE *out, PDB *pdb, ZONE *cdrs, REAL dist)
   ------------------------------------------------------------
   16.11.14   Original   By: ACRM
   03.12.14   Added output or amino acid type
*/
void RunAnalysis(FILE *out, PDB *pdb, ZONE *cdrs, REAL dist)
{
   PDB *p        = NULL,
       *nextPRes = NULL,
       *q        = NULL,
       *nextQRes = NULL;
   REAL distSq   = dist * dist,
        theDist  = 0.0;
   char theCDR[MAXLABEL];
   
   
   for(p=pdb; p!=NULL; p=nextPRes)
   {
      if(InCDR(p, cdrs, theCDR))
      {
         for(q=pdb; q!=NULL; q=nextQRes)
         {
            if(!InCDR(q, cdrs, NULL))
            {
               if((theDist=MakesContact(p, q, distSq)) >= (REAL)0.0)
               {
                  fprintf(out, "%s %s%d%s %s contacts %s%d%s %s : %.2f\n",
                          theCDR, 
                          p->chain, p->resnum, p->insert, p->resnam,
                          q->chain, q->resnum, q->insert, q->resnam,
                          theDist);
               }
            }
            nextQRes = blFindNextResidue(q);
         }
      }
      nextPRes = blFindNextResidue(p);
   }
}


/************************************************************************/
/*>REAL MakesContact(PDB *res1, PDB *res2, REAL distSq)
   ----------------------------------------------------
   16.11.14   Original   By: ACRM
*/
REAL MakesContact(PDB *res1, PDB *res2, REAL distSq)
{
   PDB *nextPRes = NULL,
       *nextQRes = NULL,
       *p, *q;
   REAL minDistSq  = (REAL)(100000.0);
   
   nextPRes = blFindNextResidue(res1);
   nextQRes = blFindNextResidue(res2);

   for(p=res1; p!=nextPRes; NEXT(p))
   {
      for(q=res2; q!=nextQRes; NEXT(q))
      {
         REAL d;
         d = DISTSQ(p, q);
         if((d > MAXBONDDISTSQ) && (d < minDistSq))
         {
            minDistSq = d;
         }
      }
   }

   if(minDistSq <= distSq)
   {
      return(sqrt(minDistSq));
   }
   
   return((REAL)-1.0);
}


/************************************************************************/
/*>BOOL InCDR(PDB *p, ZONE *cdrs, char *theCDR)
   --------------------------------------------
   16.11.14   Original   By: ACRM
   03.12.14   Sets theCDR
*/
BOOL InCDR(PDB *p, ZONE *cdrs, char *theCDR)
{
   int i;
   for(i=0; cdrs[i].start[0]; i++)
   {
      if(InPDBZoneSpec(p, cdrs[i].start, cdrs[i].stop))
      {
         if(theCDR != NULL)
            strncpy(theCDR, cdrs[i].cdr, MAXLABEL);
         return(TRUE);
      }
   }
   return(FALSE);
}


/************************************************************************/
/*>BOOL InPDBZoneSpec(PDB *p, char *start, char *stop)
   ---------------------------------------------------
   16.11.14   Original   By: ACRM
*/
BOOL InPDBZoneSpec(PDB *p, char *start, char *stop)
{
   int  startResnum,    stopResnum;
   char startChain[8],  stopChain[8],
        startInsert[8], stopInsert[8];
   
   blParseResSpec(start, startChain, &startResnum, startInsert);
   blParseResSpec(stop,  stopChain,  &stopResnum,  stopInsert);
   
   return(blInPDBZone(p, startChain, startResnum, startInsert,
                                     stopResnum,  stopInsert));
}



/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  REAL *dist, BOOL *doResList)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Returns: BOOL                Success?

   Parse the command line
   
   22.07.96 Original    By: ACRM
   29.09.05 Added uppercaseresspec param and handling of -l  By: TL
   05.11.07 Added first check that at least one parameter is on the
            command line
   04.12.14 Added doResList
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  REAL *dist, BOOL *doResList)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   *dist     = DEF_DISTCUT;

   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'h':
            return(FALSE);
            break;
         case 'r':
            *doResList = TRUE;
            break;
         case 'd':
            argc--;
            argv++;
            sscanf(argv[0], "%lf", dist);
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
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
   16.11.14   Original   By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\n");
   fprintf(stderr,"cdrcontacts V1.1 (c) 2014, Dr. Andrew C.R. Martin\n");

   fprintf(stderr,"\nUsage: cdrcontacts [-d dist] [-r] [in.pdb \
[outfile]]\n");
   fprintf(stderr,"       -d  Specify contact distance [Default: %.2f]\n",
           DEF_DISTCUT);
   fprintf(stderr,"       -r  Output the complete residue list too\n");

   fprintf(stderr,"\nFind contacts between CDRs and framework in an \
antibody\n\n");
}
