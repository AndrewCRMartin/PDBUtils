/************************************************************************/
/**

   Program:    saltbridge
   \file       saltbridge.c
   
   \version    V0.1
   \date       26.03.18   
   \brief      Identify salt bridges and, optionally, CA distance
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 2018
   \author     Dr. Andrew C. R. Martin
   \par
               Institute of Structural & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   \par
               andrew@bioinf.org.uk
               andrew.martin@ucl.ac.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified.

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
#include <math.h>
#include "bioplib/MathType.h"
#include "bioplib/macros.h"
#include "bioplib/pdb.h"
#include "bioplib/general.h"
#include "bioplib/access.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF    256
#define MAXSBATOMS   8
#define SBDISTSQ   (REAL)16.0
#define DATADIR    "DATADIR"
#define DEF_RESRADFILE "radii.dat"

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  BOOL *doHis, BOOL *doSA);
void Usage(void);
void CalculateAndDisplaySaltbridges(FILE *out, PDB *pdb, BOOL doHis,
                                    BOOL doSolv);
BOOL PrintSaltBridge(FILE *out, PDB *p, PDB *q, BOOL doHis, BOOL doSolv);
int FindSBAtoms(PDB *p, PDB **atoms, int maxsbatoms);
BOOL TestSBDistance(PDB **pAtoms, int NAtomsP, PDB **qAtoms, int NAtomsQ,
                    PDB **pAtom, PDB **qAtom, REAL *distance);
REAL CalcCADistance(PDB *p, PDB *q);
REAL SumSolv(PDB *start);



/************************************************************************/
int main(int argc, char **argv)
{
   char inFile[MAXBUFF],
        outFile[MAXBUFF];
   FILE *in  = stdin,
        *out = stdout;
   BOOL doHis = FALSE,
      doSA = FALSE,
      noEnv = FALSE;
   int  natoms;
   PDB  *pdb;
   char resradFile[MAXBUFF];

   strcpy(resradFile, DEF_RESRADFILE);
   
   if(ParseCmdLine(argc, argv, inFile, outFile, &doHis, &doSA))
   {
      if(blOpenStdFiles(inFile, outFile, &in, &out))
      {
         if((pdb=blReadPDB(in, &natoms))!=NULL)
         {
            if(doSA)
            {
               FILE   *fpRad  = NULL;

               if((fpRad = blOpenFile(resradFile, DATADIR, "r", &noEnv))==NULL)
               {
                  fprintf(stderr,"pdbsaltbridge: Unable to open residue radius file (%s)\n", resradFile);
                  if(noEnv)
                  {
                     fprintf(stderr,"               Environment variable (%s) notg set\n", DATADIR);
                  }
                  return(1);
               }
               
               blSetAtomRadii(pdb, fpRad);
               blCalcAccess(pdb, natoms, (REAL)0.0, (REAL)1.4, TRUE);
            }

            CalculateAndDisplaySaltbridges(out, pdb, doHis, doSolv);
         }
         else
         {
            fprintf(stderr,"pdbsaltbridge: Error - no atoms read from PDB \
file\n");
            return(1);
         }
      }
      else
      {
         fprintf(stderr,"pdbsaltbridge: Error - unable to open input or \
output file\n");
         return(1);
      }
   }
   else
   {
      Usage();
   }

   return(0);
}



void CalculateAndDisplaySaltbridges(FILE *out, PDB *pdb, BOOL doHis,
                                    BOOL doSolv)
{
   PDB *p, *pNextRes,
       *q, *qNextRes;
   
   for(p=pdb; p!=NULL; p=pNextRes)
   {
      pNextRes = blFindNextResidue(p);
      
      for(q=pNextRes; q!=NULL; q=qNextRes)
      {
         qNextRes = blFindNextResidue(q);
         PrintSaltBridge(out, p, q, doHis, doSolv);
      }
   }
}


BOOL PrintSaltBridge(FILE *out, PDB *p, PDB *q, BOOL doHis, BOOL doSolv)
{
   BOOL potentialSB = FALSE;
   
   if(!strncmp(p->resnam, "GLU", 3) ||
      !strncmp(p->resnam, "ASP", 3))
   {
      if(!strncmp(q->resnam, "ARG", 3) ||
         !strncmp(q->resnam, "LYS", 3) ||
         (doHis && !strncmp(q->resnam, "HIS", 3)))
      {
         potentialSB = TRUE;
      }
   }
   else if(!strncmp(p->resnam, "ARG", 3) ||
           !strncmp(p->resnam, "LYS", 3) ||
           (doHis && !strncmp(p->resnam, "HIS", 3)))
   {
      if(!strncmp(q->resnam, "GLU", 3) ||
         !strncmp(q->resnam, "ASP", 3))
      {
         potentialSB = TRUE;
      }
   }

   if(potentialSB && doSolv)
   {
      REAL totalP = (REAL)0.0,
         totalQ = (REAL)0.0;

      totalP = SumSolv(p);
      totalQ = SumSolv(q);

      if((totalP < (REAL)20.0) ||
         (totalQ < (REAL)20.0))
      {
         potentialSB = FALSE;
      }
   }

   if(potentialSB)
   {
      PDB *pAtoms[MAXSBATOMS], *qAtoms[MAXSBATOMS],
          *pAtom = NULL, *qAtom = NULL;
      REAL distance;
      int NAtomsP, NAtomsQ;

#ifdef DEBUG
      fprintf(stderr,"Testing %s %s%d%s : %s %s%d%s\n",
              p->resnam, p->chain, p->resnum, p->insert,
              q->resnam, q->chain, q->resnum, q->insert);
#endif

      NAtomsP = FindSBAtoms(p, pAtoms, MAXSBATOMS);
#ifdef DEBUG
      {
         int i;
         for(i=0; i<NAtomsP; i++)
         {
            if(pAtoms[i] != NULL)
            {
               fprintf(stderr, "%s %s%d%s %s | ",
                       pAtoms[i]->resnam,
                       pAtoms[i]->chain, pAtoms[i]->resnum, pAtoms[i]->insert,
                       pAtoms[i]->atnam);
            }
         }
         fprintf(stderr,"\n");
      }
#endif

      NAtomsQ = FindSBAtoms(q, qAtoms, MAXSBATOMS);
#ifdef DEBUG
      {
         int i;
         for(i=0; i<NAtomsQ; i++)
         {
            if(qAtoms[i] != NULL)
            {
               fprintf(stderr, "%s %s%d%s %s | ",
                       qAtoms[i]->resnam,
                       qAtoms[i]->chain, qAtoms[i]->resnum, qAtoms[i]->insert,
                       qAtoms[i]->atnam);
            }
         }
         fprintf(stderr,"\n\n");
      }
#endif

      if(TestSBDistance(pAtoms, NAtomsP, qAtoms, NAtomsQ,
                        &pAtom, &qAtom, &distance))
      {
         REAL CADistance = CalcCADistance(p, q);
         
         fprintf(out,
                 "SB: %s %s%d%s %s | %s %s%d%s %s Dist: %.3f CA: %.3f\n",
                 pAtom->resnam,
                 pAtom->chain, pAtom->resnum, pAtom->insert, pAtom->atnam,
                 qAtom->resnam,
                 qAtom->chain, qAtom->resnum, qAtom->insert, qAtom->atnam,
                 distance, CADistance);
         return(TRUE);
      }
   }
   
   return(FALSE);
}


REAL SumSolv(PDB *start)
{
   PDB *p, *stop;
   REAL sum = (REAL)0.0;
   
   stop = blFindNextResidue(start);
   for(p=start; p!=stop; NEXT(p))
   {
      sum += p->access;
   }
   return(sum);
   
}

REAL CalcCADistance(PDB *p, PDB *q)
{
   PDB *pCA = NULL,
      *qCA = NULL;
   pCA = blFindAtomInRes(p, "CA  ");
   qCA = blFindAtomInRes(q, "CA  ");

   if((pCA!=NULL) && (qCA != NULL))
   {
      return(DIST(pCA, qCA));
   }
   
   return((REAL)0.0);
}



/* Note return values must be < MAXSBATOMS */
int FindSBAtoms(PDB *p, PDB **atoms, int maxsbatoms)
{
   int i;
   for(i=0; i<maxsbatoms; i++)
   {
      atoms[i] = NULL;
   }
   if(!strncmp(p->resnam, "ASP", 3))
   {
      atoms[0] = blFindAtomInRes(p, "OD1");
      atoms[1] = blFindAtomInRes(p, "OD2");
      return(2);
   }
   else if(!strncmp(p->resnam, "GLU", 3))
   {
      atoms[0] = blFindAtomInRes(p, "OE1");
      atoms[1] = blFindAtomInRes(p, "OE2");
      return(2);
   }
   else if(!strncmp(p->resnam, "ARG", 3))
   {
      atoms[0] = blFindAtomInRes(p, "NE");
      atoms[1] = blFindAtomInRes(p, "NH1");
      atoms[2] = blFindAtomInRes(p, "NH2");
      return(3);
   }
   else if(!strncmp(p->resnam, "LYS", 3))
   {
      atoms[0] = blFindAtomInRes(p, "NZ");
      return(1);
   }
   else if(!strncmp(p->resnam, "HIS", 3))
   {
      atoms[0] = blFindAtomInRes(p, "ND1");
      atoms[1] = blFindAtomInRes(p, "NE2");
      return(2);
   }
   return(0);
}


BOOL TestSBDistance(PDB **pAtoms, int NAtomsP, PDB **qAtoms, int NAtomsQ,
                    PDB **pAtom, PDB **qAtom, REAL *distance)
{
   int CountP, CountQ;
   REAL bestDistSq = SBDISTSQ;
   BOOL found = FALSE;

   /* Value of zero indicates no distance found */
   *distance = (REAL)0.0;

   /* Step through pairs of atoms */
   for(CountP=0; CountP<NAtomsP; CountP++)
   {
      for(CountQ=0; CountQ<NAtomsQ; CountQ++)
      {
         REAL distSq = (REAL)0.0;
         if((pAtoms[CountP] != NULL) &&
            (qAtoms[CountQ] != NULL))
         {
            distSq = DISTSQ(pAtoms[CountP], qAtoms[CountQ]);

            /* Keep the closest pair */
            if(distSq <= bestDistSq)
            {
               *pAtom = pAtoms[CountP];
               *qAtom = qAtoms[CountQ];
               bestDistSq = distSq;
               found = TRUE;
            }
         }
      }
   }

   /* If we found a pair within required distance, calculate the actual 
      distance
   */
   if(found)
   {
      *distance = sqrt(bestDistSq);
   }
   
   return(found);
}




/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                     BOOL *CATorsions, BOOL *terse, BOOL *Radians, 
                     BOOL *oldStyle)
   ---------------------------------------------------------------------
*//**

   \param[in]     argc         Argument count
   \param[in]     **argv       Argument array
   \param[out]    *infile      Input file (or blank string)
   \param[out]    *outfile     Output file (or blank string)
   \param[out]    *CATorsions  Do CA pseudo-torsions
   \param[out]    *terse       Terse (1-letter AA code) output
   \param[out]    *Radians     Output radians rather than degrees
   \param[out]    *oldStyle    Old style output
   \return                     Success?

   Parse the command line
   
-  05.02.96 Original    By: ACRM
-  27.02.14 V2.0
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  BOOL *doHis, BOOL *doSA)
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
         case 'H':
            *doHis = TRUE;
            break;
         case 's':
            *doSA = TRUE;
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are <= 2 arguments left                    */
         if(argc > 2)
            return(FALSE);
         
         /* Copy the first to infile                                    */
         if(argc)
         {
            strcpy(infile, argv[0]);
            argc--;
         }
         
         /* Copy the second to outfile                                  */
         if(argc)
         {
            strcpy(infile, argv[0]);
            argc--;
         }
         
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}

void Usage(void)
{

}

