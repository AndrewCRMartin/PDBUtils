/************************************************************************/
/**

   Program:    protrusion
   \file       protrusion.c
   
   \version    V1.0
   \date       12.01.21   
   \brief         
   
   \copyright  (c) UCL / Prof. Andrew C. R. Martin 2021
   \author     Prof. Andrew C. R. Martin
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
#include <math.h>
#include "bioplib/macros.h"
#include "bioplib/pdb.h"
#include "bioplib/MathType.h"
#include "bioplib/MathUtil.h"

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
PDB *RunAnalysis(PDB *pdb, char *startres, char *stopres,
                 REAL *protrusion);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  char *startres, char *stopres);
void Usage(void);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**
   Main program

- 12.01.21 Original   By: ACRM
**/
int main(int argc, char **argv)
{
   FILE *in     = stdin,
        *out    = stdout;
   int  natoms;
   REAL protrusion = 0.0;
   PDB  *pdb, *protrudingRes;
   char infile[MAXBUFF],
        outfile[MAXBUFF],
        startres[MAXBUFF],
        stopres[MAXBUFF];
   
   if(!ParseCmdLine(argc, argv, infile, outfile,
                    startres, stopres))
   {
      Usage();
      return(0);
   }

   if(!blOpenStdFiles(infile, outfile, &in, &out))
   {
      fprintf(stderr, "Error (protrusion): Unable to open input or output \
file\n");
      return(1);
   }

   if((pdb = blReadPDB(in, &natoms))==NULL)
   {
      fprintf(stderr, "Error (protrusion): No atoms read from PDB \
file, %s\n", infile);
      return(1);
   }

   /* Reduce to CA atoms only                                           */
   pdb = blSelectCaPDB(pdb);

   if((protrudingRes=RunAnalysis(pdb, startres, stopres, &protrusion))
      != NULL)
   {
      fprintf(out, "%s%d%s %s Protrusion: %.3f\n",
              protrudingRes->chain,
              protrudingRes->resnum,
              protrudingRes->insert,
              protrudingRes->resnam,
              protrusion);
   }
   
   return(0);
}

/************************************************************************/
/*>PDB *RunAnalysis(PDB *pdb, char *startres, char *stopres,
                    REAL *protrusion)
   --------------------------------------------------------------------
*//**
   \param[in]   PDB    *pdb             PDB linked list
   \param[in]   char   *startres        start residue specification
   \param[in]   char   *stopres         stop residue specification
   \param[out]  REAL   *protrusion      The protrusion of the most
                                        protruding residue
   \return      PDB*                    The PDB record of the most
                                        protruding residue

   Calculates the protrusion of the most protrtuding residue from the
   line between the two specified ends

-  12.01.21 Original   By: ACRM
**/
PDB *RunAnalysis(PDB *pdb, char *startres, char *stopres,
                 REAL *protrusion)
{
   PDB   *p, *start, *stop,
         *residue=NULL;
   VEC3F end1, end2;
   REAL  maxDist = 0.0;

   /* Find the start and stop of the zone of interest                   */
   if((start = blFindResidueSpec(pdb, startres))==NULL)
      return(NULL);
   if((stop  = blFindResidueSpec(pdb, stopres))==NULL)
      return(NULL);

   /* Fill in the vectors                                               */
   end1.x = start->x;
   end1.y = start->y;
   end1.z = start->z;
   
   end2.x = stop->x;
   end2.y = stop->y;
   end2.z = stop->z;
   
   /* Step through the PDB linked list between start and stop and 
      measure the distance from the line, recording the longest
   */
   for(p=start->next; p!=stop; NEXT(p))
   {
      VEC3F point;
      REAL  dist;

      point.x = p->x;
      point.y = p->y;
      point.z = p->z;
      
      dist = blDistPtLine(point, end1, end2);
      if(dist > maxDist)
      {
         maxDist = dist;
         residue = p;
      }
   }

   *protrusion = maxDist;
   return(residue);
}
      

/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                     char *startres, char *stopres)
   ----------------------------------------------------------------------
*//**
   \param[in]   int    argc              Argument count
   \param[in]   char   **argv            Argument array
   \param[out]  char   *infile           Input filename (or blank string)
   \param[out]  char   *outfile          Output filename (or blank string)
   \param[out]  char   *startres         First residue of zone
   \param[out]  char   *stopres          Last residue of zone

   \return      BOOL                     Success

   Parse the command line

   17.07.14 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  char *startres, char *stopres)
{
   argc--;
   argv++;
   
   infile[0] = outfile[0] = startres[0] = stopres[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'h':
            return(FALSE);
         default:
            return(FALSE);
         }
      }
      else
      {
         /* Check that there are between 2 and 4 arguments left         */
         if((argc < 2) || (argc > 4))
            return(FALSE);
         
         /* Copy the first and secont to startres/stopres               */
         strcpy(startres, argv[0]);
         argc--;
         argv++;
         strcpy(stopres,  argv[0]);
         argc--;
         argv++;
         
         /* If there's another, copy it to infile                      */
         if(argc)
         {
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
*//**
   Prints a usage message

-   12.01.21 Original   By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\nprotrusion V1.1 (c) 2021 UCL, Prof. Andrew C.R. \
Martin\n");

   fprintf(stderr,"\nUsage: protrusion startres lastres \
[in.pdb [out.txt]]\n");

   fprintf(stderr,"\nTakes a zone of residues and finds the distance of \
the intervening\n");
   fprintf(stderr,"residue that protrudes most from the line between \
the two specified\n");
   fprintf(stderr,"residues. Uses only the C-alpha atoms.\n\n");
}

