/*************************************************************************

   Program:    splitloop
   File:       splitloop.c
   
   Version:    V1.2
   Date:       24.10.94
   Function:   Take a PDB file of a section of protein purporting to be
               a loop. Divides it up into real loops (i.e. sections
               obeying a set of rules defined below) and outputs these
               separated by TER cards.
   
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

   Regions of a PDB file which are 'coil' may involve a number of actual
   loops. This program devides such sections (supplied as an input PDB
   file) into (possibly overlapping) sections which obey the following
   rules.

   1. The endpoints are separated by less than SPANFRAC % of the total 
      possible span of the loop.
   2. The loop does not cross a vector between its ends (i.e. it lies
      totally on one side of the vector).
   3. The CofG of CAs is at least MINCG Angstroms from the endpoint
      vector.

   N.B. After this procedure, the loop(s) will have been moved in space.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  20.10.94 Original    By: ACRM
   V1.1  21.10.94 Modified default span to 65%. Added -l option
   V1.2  24.10.94 Modified to return wider loops rather than look for
                  the nearer CA

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/matrix.h"
#include "bioplib/general.h"
#include "bioplib/angle.h"

/************************************************************************/
/* Defines and macros
*/
#define CACADIST 3.8
#define SPANFRAC 0.65
#define MINCG    2.0

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL AnalyseLoop(FILE *out, PDB *pdb, REAL SpanFrac, REAL MinCG, 
                 int MinLen);
BOOL EndPointsOK(PDB *pdb, REAL SpanFrac);
void FindCAlphas(PDB *pdb, VEC3F *nter, VEC3F *cter, VEC3F *secres);
void OrientLoop(PDB *pdb);
void RotateToXZ(PDB *pdb, VEC3F *Cter, VEC3F *CofG);
void RotateToX(PDB *pdb, VEC3F *Cter, VEC3F *CofG);
void RotateLoopToXY(PDB *pdb, VEC3F *Cter, VEC3F *CofG);
BOOL IsDoubleLoop(PDB *pdb, PDB **end, PDB **start);
PDB *BuildFirstPDB(PDB *pdb, PDB *end);
PDB *BuildSecondPDB(PDB *pdb, PDB *start);
BOOL LoopOK(PDB *pdb, REAL MinCG, int MinLen);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,    
                  REAL *SpanFrac, REAL *MinCG, int *MinLen);
void Usage(void);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for loop analysis.

   20.10.94 Original    By: ACRM
   21.10.94 Added MinLen
*/
int main(int argc, char **argv)
{
   FILE *in      = stdin,
        *out     = stdout;
   PDB  *pdb;
   int  natom,
        MinLen   = 0;
   REAL SpanFrac = SPANFRAC,
        MinCG    = MINCG;
   char infile[160],
        outfile[160];

   if(ParseCmdLine(argc, argv, infile, outfile, &SpanFrac, &MinCG, 
                   &MinLen))
   {
      if(OpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdb=ReadPDB(in, &natom)) != NULL)
         {
            if(!AnalyseLoop(out, pdb, SpanFrac, MinCG, MinLen))
               return(1);
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
/*>BOOL AnalyseLoop(FILE *out, PDB *pdb, REAL SpanFrac, REAL MinCG, 
                    int MinLen)
   ----------------------------------------------------------------
   Recursive routine to perform the analysis of a loop. Writes out a
   PDB file to file 'out' if all OK. Otherwise recurses.

   Returns FALSE only if there is a memory error.

   20.10.94 Original    By: ACRM
   21.10.94 Added MinLen
*/
BOOL AnalyseLoop(FILE *out, PDB *pdb, REAL SpanFrac, REAL MinCG, 
                 int MinLen)
{
   PDB *pdb1,
       *pdb2,
       *end,
       *start;

   if(EndPointsOK(pdb, SpanFrac))
   {
      OrientLoop(pdb);

      if(IsDoubleLoop(pdb, &end, &start))
      {
         /* Split the loop into the two parts                           */
         if((pdb1 = BuildFirstPDB(pdb, end))==NULL)
            return(FALSE);
         if((pdb2 = BuildSecondPDB(pdb, start))==NULL)
            return(FALSE);

         /* Free the memory used by the input PDB linked list           */
         FREELIST(pdb, PDB);

         /* Recurse to analyse each of these loop sections              */
         if(!AnalyseLoop(out, pdb1, SpanFrac, MinCG, MinLen))
            return(FALSE);
         if(!AnalyseLoop(out, pdb2, SpanFrac, MinCG, MinLen))
            return(FALSE);
      }
      else         /* It's a single loop                                */
      {
         if(LoopOK(pdb, MinCG, MinLen))
         {
            WritePDB(out, pdb);
         }
      }
   }
   
   return(TRUE);
}

/************************************************************************/
/*>BOOL EndPointsOK(PDB *pdb, REAL SpanFrac)
   -----------------------------------------
   Checks whether the terminal CAs are less than SpanFrac * max possible
   separation. i.e. the segment isn't extended

   20.10.94 Original    By: ACRM
   21.10.94 Fixed MaxSpan calc to NRes-1
*/
BOOL EndPointsOK(PDB *pdb, REAL SpanFrac)
{
   PDB  *FirstCA = NULL, 
        *LastCA  = NULL,
        *p;
   int  NRes     = 0;
   REAL MaxSpan,
        Span;

   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(!strncmp(p->atnam, "CA  ", 4))
      {
         if(FirstCA == NULL)
            FirstCA = p;
         LastCA = p;
         NRes++;
      }
   }

   /* Calculate theoretical maximum span and the actual span            */
   MaxSpan = (NRes - 1) * CACADIST;
   Span = DIST(FirstCA, LastCA);

   if((Span/MaxSpan) <= SpanFrac)
      return(TRUE);
   
   return(FALSE);
}

/************************************************************************/
/*>void FindCAlphas(PDB *pdb, VEC3F *nter, VEC3F *cter, VEC3F *secres)
   -------------------------------------------------------------------
   Finds coordinates for the first, third and last CAs

   20.10.94 Original    By: ACRM
   24.10.94 Uses third rather than second CA if there is one
*/
void FindCAlphas(PDB *pdb, VEC3F *nter, VEC3F *cter, VEC3F *secres)
{
   PDB  *FirstCA  = NULL, 
        *SecondCA = NULL,
        *LastCA   = NULL,
        *p;
   int  NRes      = 0;

   /* Find pointers to the three CAs                                    */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(!strncmp(p->atnam, "CA  ", 4))
      {
         if(FirstCA == NULL)
            FirstCA = p;

         LastCA = p;

         if(++NRes < 4)      /* i.e. is done for CAs 1, 2, 3            */
            SecondCA = p;
      }
   }

   /* Set the 3 vectors                                                 */
   nter->x   = FirstCA->x;
   nter->y   = FirstCA->y;
   nter->z   = FirstCA->z;

   cter->x   = LastCA->x;
   cter->y   = LastCA->y;
   cter->z   = LastCA->z;

   secres->x = SecondCA->x;
   secres->y = SecondCA->y;
   secres->z = SecondCA->z;
}


/************************************************************************/
/*>void OrientLoop(PDB *pdb)
   -------------------------
   Orients a loop such that the terminal CAs are along the x-axis and
   the CA or residue 2 is on the xy-plane

   20.10.94 Original    By: ACRM
*/
void OrientLoop(PDB *pdb)
{
   VEC3F SecRes,
         Nter,
         Cter,
         TempVec;

   FindCAlphas(pdb, &Nter, &Cter, &SecRes);
   
   /* Move the PDB so Nter is at the origin                             */
   TempVec.x = -Nter.x;
   TempVec.y = -Nter.y;
   TempVec.z = -Nter.z;
   
   TranslatePDB(pdb, TempVec);

   /* Refind Coords of terminii since Cter has moved                    */
   FindCAlphas(pdb, &Nter, &Cter, &SecRes);
   
   /* Rotate the Cter onto the XZ plane                                 */
   RotateToXZ(pdb, &Cter, &SecRes);
   
   /* Rotate the Cter onto the X axis                                   */
   RotateToX(pdb, &Cter, &SecRes);
   
   /* Now rotate about the X axis such that CofG is on the XY plane     */
   RotateLoopToXY(pdb, &Cter, &SecRes);
}

/************************************************************************/
/*>void RotateToXZ(PDB *pdb, VEC3F *Cter, VEC3F *CofG)
   ---------------------------------------------------
   I/O:     PDB    *pdb           PDB linked list
            VEC3F  *Cter          Cter CA coordinates

   Rotate loop such that Cter CA is on the XZ plane

   25.07.94 Original    By: ACRM
   29.07.94 Corrected Rotate calls to rotate both Cter and CofG
*/
void RotateToXZ(PDB *pdb, VEC3F *Cter, VEC3F *CofG)
{
   REAL  ang;
   REAL  matrix[3][3];
   VEC3F OutVec;
   
   ang = TrueAngle(Cter->y, Cter->x);
   
   CreateRotMat('z', -ang, matrix);
   
   ApplyMatrixPDB(pdb, matrix);
   MatMult3_33(*Cter,matrix,&OutVec);
   *Cter = OutVec;
   MatMult3_33(*CofG,matrix,&OutVec);
   *CofG = OutVec;
}

/************************************************************************/
/*>void RotateToX(PDB *pdb, VEC3F *Cter, VEC3F *CofG)
   --------------------------------------------------
   I/O:     PDB    *pdb           PDB linked list
            VEC3F  *Cter          Cter CA coordinates

   Having called RotateToXZ(), rotate loop such that Cter CA is on the
   X axis

   25.07.94 Original    By: ACRM
   29.07.94 Corrected Rotate calls to rotate both Cter and CofG
*/
void RotateToX(PDB *pdb, VEC3F *Cter, VEC3F *CofG)
{
   REAL  ang;
   REAL  matrix[3][3];
   VEC3F OutVec;
   
   ang = TrueAngle(Cter->z, Cter->x);
   
   CreateRotMat('y', ang, matrix);
   
   ApplyMatrixPDB(pdb, matrix);
   MatMult3_33(*Cter,matrix,&OutVec);
   *Cter = OutVec;
   MatMult3_33(*CofG,matrix,&OutVec);
   *CofG = OutVec;
}

/************************************************************************/
/*>void RotateLoopToXY(PDB *pdb, VEC3F *Cter, VEC3F *CofG)
   -------------------------------------------------------
   I/O:     PDB    *pdb           PDB linked list
            VEC3F  *CofG          CofG coordinates

   Having orientated loop along the x-axis, rotate such that CofG is on
   the XY plane

   25.07.94 Original    By: ACRM
   29.07.94 Corrected Rotate calls to rotate both Cter and CofG
*/
void RotateLoopToXY(PDB *pdb, VEC3F *Cter, VEC3F *CofG)
{
   REAL  ang;
   REAL  matrix[3][3];
   VEC3F OutVec;
   
   ang = TrueAngle(CofG->z, CofG->y);
   
   CreateRotMat('x', -ang, matrix);
   
   ApplyMatrixPDB(pdb, matrix);
   MatMult3_33(*Cter,matrix,&OutVec);
   *Cter = OutVec;
   MatMult3_33(*CofG,matrix,&OutVec);
   *CofG = OutVec;
}


/************************************************************************/
/*>BOOL IsDoubleLoop(PDB *pdb, PDB **end, PDB **start)
   ---------------------------------------------------
   Once the loop has been correctly oriented along the x-axis, any
   negative C-alpha y-coordinates indicate a double loop. The routine
   returns TRUE if this is the case and outputs PDB pointers for the 
   last CA in the first sub-loop and the first CA in the second sub-loop.

   20.10.94 Original    By; ACRM
   24.10.94 Modified to return wider loop rather than to nearer CA
            Makes a check that the split loops don't actually encompass
            the whole loop (and if so, returns FALSE)
*/
BOOL IsDoubleLoop(PDB *pdb, PDB **end, PDB **start)
{
   PDB *p,
       *FirstCA = NULL,
       *PrevCA  = NULL,
       *LastCA  = NULL;

   /* Find the last CA                                                  */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(!strncmp(p->atnam,"CA  ",4))
      {
         LastCA = p;
      }
   }

#ifdef DEBUG
   fprintf(stderr,"Testing Loop from %d to %d\n",
           pdb->resnum,
           LastCA->resnum);
#endif   
   
   /* Run through the CAs in the linked list to see if we cross the
      vector between the termini
   */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(!strncmp(p->atnam,"CA  ",4))
      {
         /* Store the first CA                                          */
         if(FirstCA == NULL)
            FirstCA = p;
         
         /* If y is negative, we have crossed the vector                */
         if(p->y < (REAL)(-0.5))
         {
            /* See if this CA or the previous CA is closer to the N-ter
               CA and to the C-ter CA
            */
            if(PrevCA == NULL)
            {
               *end = *start = p;
            }
            else
            {
               /* Set sub-loops to wider possible spans                 */
               *end   = p;
               *start = PrevCA;

               /* Test that these don't actually cover the whole loop   */
               if(*end == LastCA ||
                  *start == FirstCA)
                  return(FALSE);
               
#ifdef DEBUG
               fprintf(stderr,"Subloops: %d-%d and %d-%d\n",
                       pdb->resnum, (*end)->resnum,
                       (*start)->resnum, LastCA->resnum);
#endif

               
               /* First the end of the section from the N-terminus      */
/*
//               if(DISTSQ(FirstCA, p) < DISTSQ(FirstCA, PrevCA))
//                  *end = p;
//               else
//                  *end = PrevCA;
*/
               /* Now the start of the section ending at the C-terminus */
/*
//               if(DISTSQ(LastCA, p) < DISTSQ(LastCA, PrevCA))
//                  *start = p;
//               else
//                  *start = PrevCA;
*/
            }

            /* Indicate that this was a double loop                     */
            return(TRUE);
         }
         PrevCA = p;
      }
   }

   /* Was not a double loop                                             */
   return(FALSE);
}

/************************************************************************/
/*>PDB *BuildFirstPDB(PDB *pdb, PDB *end)
   --------------------------------------
   Build and return a new linked list to run from the start of the 
   original list to the end of the residue containing end

   20.10.94 Original    By: ACRM
*/
PDB *BuildFirstPDB(PDB *pdb, PDB *end)
{
   PDB  *p, *q,
        *retpdb   = NULL;
   BOOL GotEnd    = FALSE;
   char PrevChain = '-', 
        PrevIns   = '-';
   int  PrevRes   = -1000;

   for(p=pdb; p!=NULL; NEXT(p))
   {
      /* Set flag if we've got to the atom whose residue we need to end
         at
      */
      if(p==end) GotEnd = TRUE;

      /* See if we go beyond the residue containing the atom at which
         we need to end
      */
      if(GotEnd && ((p->resnum    != PrevRes) ||
                    (p->insert[0] != PrevIns) ||
                    (p->chain[0]  != PrevChain)))
      {
         return(retpdb);
      }

      /* Allocate space in new linked list                              */
      if(retpdb==NULL)
      {
         INIT(retpdb,PDB);
         q = retpdb;
      }
      else
      {
         ALLOCNEXT(q, PDB);
      }
      
      /* Check allocation                                               */
      if(q==NULL)
         return(NULL);

      /* Copy item into new linked list                                 */
      CopyPDB(q, p);

      /* Update record of previous entry in linked list                 */
      PrevRes   = p->resnum;
      PrevIns   = p->insert[0];
      PrevChain = p->chain[0];
   }

   return(retpdb);
}

/************************************************************************/
/*>PDB *BuildSecondPDB(PDB *pdb, PDB *start)
   -----------------------------------------
   Builds a PDB linked list starting at the first atom of the residue
   containing start

   20.10.94 Original    By: ACRM
*/
PDB *BuildSecondPDB(PDB *pdb, PDB *start)
{
   PDB *p, *q,
       *thisres,
       *nextres,
       *retpdb = NULL;

   /* For each residue in the PDB linked list                           */
   for(thisres = pdb; thisres != NULL; thisres = nextres)
   {
      /* Find the start of the next residue                             */
      nextres = FindEndPDB(thisres);

      /* See if our start atom is in this residue                       */
      for(p=thisres; p!=nextres; NEXT(p))
      {
         if(p==start)    /* Found it                                    */
         {
            /* Just copy the linked list from thisres on to the end     */
            for(p=thisres; p!=NULL; NEXT(p))
            {
               /* Allocate space in new linked list                     */
               if(retpdb==NULL)
               {
                  INIT(retpdb,PDB);
                  q = retpdb;
               }
               else
               {
                  ALLOCNEXT(q, PDB);
               }
               
               /* Check allocation                                      */
               if(q==NULL)
                  return(NULL);
               
               /* Copy item into new linked list                        */
               CopyPDB(q, p);
            }

            /* And return our new linked list                           */
            return(retpdb);
         }
      }
   }
   
   return(NULL);
}

/************************************************************************/
/*>BOOL LoopOK(PDB *pdb, REAL MinCG, int MinLen)
   ---------------------------------------------
   Checks that the loop is at least MinLen residues long.

   WILL CHECK THE DISTANCE BETWEEN THE CofG OF THE CAs AND THE VECTOR 
   BETWEEN THE TERMINAL CAs.

   21.10.94 Original    By: ACRM
*/
BOOL LoopOK(PDB *pdb, REAL MinCG, int MinLen)
{
   PDB *p;
   int NRes = 0;
   
   /* Count the CAs in the loop                                         */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(!strncmp(p->atnam, "CA  ", 4))
         NRes++;
   }

   if(NRes >= MinLen)
   {
      /* SHOULD NOW CHECK THE CofG TO VECTOR DISTANCE...                */

      return(TRUE);
   }
   
   return(FALSE);
}

/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                     REAL *SpanFrac, REAL *MinCG, int *MinLen)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
            REAL   *SpanFrac    Max allowed terminal separation
            REAL   *MinCG       Min dist of CofG from terminal vector
            int    *MinLen      Minumum length of a loop (default: 0)
   Returns: BOOL                Success?

   Parse the command line
   
   20.10.94 Original    By: ACRM
   21.10.94 Added MinLen
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  REAL *SpanFrac, REAL *MinCG, int *MinLen)
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
         case 'w':
            argc--;
            argv++;
            sscanf(argv[0],"%lf",SpanFrac);
            *SpanFrac /= (REAL)100.0;
            break;
         case 'c':
            argc--;
            argv++;
            sscanf(argv[0],"%lf",MinCG);
            break;
         case 'l':
            argc--;
            argv++;
            sscanf(argv[0],"%d",MinLen);
            break;
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

   20.10.94 Original    By: ACRM
   21.10.94 V1.1
   24.10.94 V1.2
*/
void Usage(void)
{
   fprintf(stderr,"\nSplitLoop V1.2 (c) 1994, Andrew C.R. Martin, \
UCL\n\n");
   fprintf(stderr,"Usage: splitloop [-w width] [-c CofG] [-l minlen] \
[<in.pdb>] [<out.pdb>]\n");
   fprintf(stderr,"       -w Specify percentage of max terminal \
separation allowed\n");
   fprintf(stderr,"       -c Specify minimum distance of CofG from \
terminal vector.\n");
   fprintf(stderr,"       -l Specify minimum length of a loop \
(default: 0).\n\n");
   fprintf(stderr,"Split an input PDB file supposedly containing a loop \
into true loops.\n");
   fprintf(stderr,"The endpoints must be separated by less than width \
percent of the\n");
   fprintf(stderr,"maximum possible separation. The loop will not cross \
the vector between\n");
   fprintf(stderr,"its terminal C-alphas and the CofG of the loop \
C-alphas will be at least\n");
   fprintf(stderr,"CofG Angstroms from the vector.\n\n");
}


#ifdef JUNK_CODE

/* This part was an attempt to walk to the nearest CA across the loop   
   It should be included in IsDoubleLoop() instead of just setting to
   the wider possible span

   The routine should also return FALSE if *end == LastCA
*/
               /* Start off assuming end is this CA and start is the
                  previous CA
               */
               *end   = p;
               *start = PrevCA;

               /* Find the mean of the CA positions                     */
               MeanCA.x = (p->x + PrevCA->x) / (REAL)2.0;
               MeanCA.y = (p->y + PrevCA->y) / (REAL)2.0;
               MeanCA.z = (p->z + PrevCA->z) / (REAL)2.0;

               PrevCA = p;

               /* Find the distance to the N-ter                        */
               PrevDistSq = DISTSQ(FirstCA, (&MeanCA));
               
               /* Step though the PDB linked list from here on          */
               for( ; p!=NULL; NEXT(p))
               {
                  if(!strncmp(p->atnam,"CA  ",4))
                  {
                     *end = PrevCA;

                     /* Find the mean of the CA positions               */
                     MeanCA.x = (p->x + PrevCA->x) / (REAL)2.0;
                     MeanCA.y = (p->y + PrevCA->y) / (REAL)2.0;
                     MeanCA.z = (p->z + PrevCA->z) / (REAL)2.0;

                     /* Find the distance to the N-ter                  */
                     DistSq = DISTSQ(FirstCA, (&MeanCA));

                     /* If the distance has increased, we end           */
                     if(DistSq > PrevDistSq)
                     {
                        break;
                     }
                     else     /* we continue                            */
                     {
                        PrevCA     = p;
                        PrevDistSq = DistSq;
                     }
                  }
               }
               

#endif
