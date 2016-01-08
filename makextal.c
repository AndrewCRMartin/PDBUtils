/*************************************************************************

   Program:    makextal
   File:       makextal.c
   
   Version:    V1.0
   Date:       12.10.95
   Function:   Make a crystal lattice from a PDB file
   
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
#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/general.h"
#include "bioplib/fsscanf.h"

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
void WriteLatticeStructures(FILE *out, PDB *pdb, 
                            VEC3F xtrans, VEC3F ytrans, VEC3F ztrans);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  BOOL *UserParams, VEC3F *UnitCell, VEC3F *CellAngles);
void Usage(void);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for creating a crystal lattice from a PDB file.

   11.10.95 Original   By: ACRM
*/
int main(int argc, char **argv)
{
   PDB   *pdb;
   int   natom, 
         i, j;
   VEC3F UnitCell,
         CellAngles,
         xtrans, ytrans, ztrans;
   FILE  *in  = stdin,
         *out = stdout;
   char  infile[MAXBUFF],
         outfile[MAXBUFF],
         spacegroup[MAXBUFF];
   BOOL  UserParams = FALSE;
   REAL  OrigMatrix[3][4],
         ScaleMatrix[3][4];


   /* Parse the command line                                            */
   if(ParseCmdLine(argc, argv, infile, outfile, &UserParams, &UnitCell,
                   &CellAngles))
   {
      /* Open input and output files                                    */
      if(OpenStdFiles(infile, outfile, &in, &out))
      {
         /* If unit cell parameters have not been given by the user     */
         if(!UserParams)
         {
            /* Try to get the parameters out of the PDB file            */
            if(!GetCrystPDB(in, &UnitCell, &CellAngles, spacegroup,
                            OrigMatrix, ScaleMatrix))
            {
               /* If we failed, print an error message and die          */
               fprintf(stderr,"makextal: Unit cell parameters missing \
from PDB file and not\nspecified on command line.\n");
               return(1);
            }
         }
         else
         {
            /* Crystal parameters have come from the user; set defaults
               for the Origin and Scale
            */
            for(i=0; i<3; i++)
            {
               for(j=0; j<4; j++)
               {
                  OrigMatrix[i][j]  = (REAL)0.0;
                  ScaleMatrix[i][j] = (REAL)0.0;
               }
               OrigMatrix[i][i] = (REAL)1.0;
            }
            strcpy(spacegroup,"P");
         }
         
         /* Read in the PDB file (all atoms, we're not fussy!)          */
         pdb = ReadPDBAll(in, &natom);

         /* Calculate the translation vectors to build the lattice      */
         CalcCellTrans(UnitCell, CellAngles, &xtrans, &ytrans, &ztrans);

         /* Write the CRYST record to the output file                   */
         WriteCrystPDB(out, UnitCell, CellAngles, spacegroup,
                       OrigMatrix, ScaleMatrix);
         
         /* Write the base structure out                                */
         WritePDB(out,pdb);
         
         /* Write the translated structures                             */
         WriteLatticeStructures(out, pdb, xtrans, ytrans, ztrans);
      }
   }
   else
   {
      Usage();
   }

   return(0);
}

   

/************************************************************************/
/*>void WriteLatticeStructures(FILE *out, PDB *pdb, VEC3F xtrans, 
                               VEC3F ytrans, VEC3F ztrans)
   --------------------------------------------------------------
   Input:   FILE    *out       Output file
            PDB     *pdb       PDB linked list
            VEC3F   xtrans     X translation vector
            VEC3F   ytrans     Y translation vector
            VEC3F   ztrans     Z translation vector

   Write the set of 6 surrounding structures to make the lattice.

   12.10.95 Original    By: ACRM
*/
void WriteLatticeStructures(FILE *out, PDB *pdb, 
                            VEC3F xtrans, VEC3F ytrans, VEC3F ztrans)
{
   PDB *pdb2;

   /* Make a working copy of the PDB linked list                        */
   pdb2 = DupePDB(pdb);

   /* Move along X and write out                                        */
   TranslatePDB(pdb2, xtrans);
   WritePDB(out,pdb2);

   /* Reset, move along Y and write out                                 */
   CopyPDBCoords(pdb2, pdb);
   TranslatePDB(pdb2, ytrans);
   WritePDB(out,pdb2);
   
   /* Reset, move along Z and write out                                 */
   CopyPDBCoords(pdb2, pdb);
   TranslatePDB(pdb2, ztrans);
   WritePDB(out,pdb2);

   /* Negate all the vectors                                            */
   xtrans.x *= (-1);
   xtrans.y *= (-1);
   xtrans.z *= (-1);
   ytrans.x *= (-1);
   ytrans.y *= (-1);
   ytrans.z *= (-1);
   ztrans.x *= (-1);
   ztrans.y *= (-1);
   ztrans.z *= (-1);

   /* Reset, move along X and write out                                 */
   CopyPDBCoords(pdb2, pdb);
   TranslatePDB(pdb2, xtrans);
   WritePDB(out,pdb2);

   /* Reset, move along Y and write out                                 */
   CopyPDBCoords(pdb2, pdb);
   TranslatePDB(pdb2, ytrans);
   WritePDB(out,pdb2);
   
   /* Reset, move along Z and write out                                 */
   CopyPDBCoords(pdb2, pdb);
   TranslatePDB(pdb2, ztrans);
   WritePDB(out,pdb2);

   FREELIST(pdb2, PDB);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                     BOOL *UserParams, VEC3F *UnitCell, VEC3F *CellAngles)
   -----------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Arguments
   Output:  char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
            BOOL   *UserParams  User has specified cell params
            VEC3F  *UnitCell    Unit cell dimensions
            VEC3F  *CellAngles  Unit cell angles
   Returns: BOOL                Success?

   Parse the command line

   12.10.95 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  BOOL *UserParams, VEC3F *UnitCell, VEC3F *CellAngles)
{
   int  i;
   REAL temp;
   
   argc--;
   argv++;

   *UserParams = FALSE;  
   infile[0]   = outfile[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'c':
            *UserParams = TRUE;
            for(i=0; i<6; i++)
            {
               argc--;
               argv++;
               if(argc <= 0)
                  return(FALSE);
               if(!sscanf(argv[0],"%lf",&temp))
                  return(FALSE);
               switch(i)
               {
               case 0:
                  UnitCell->x = temp;
                  break;
               case 1:
                  UnitCell->y = temp;
                  break;
               case 2:
                  UnitCell->z = temp;
                  break;
               case 3:
                  CellAngles->x = PI*temp/180.0;
                  break;
               case 4:
                  CellAngles->y = PI*temp/180.0;
                  break;
               case 5:
                  CellAngles->z = PI*temp/180.0;
                  break;
               }
            }      
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

   12.10.95 Original    By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\n\nmakextal V1.0 (c) 1995 Dr. Andrew C.R. Martin, \
UCL.\n");
   fprintf(stderr,"              Freely distributable if no profit is \
made in so doing.\n");

   fprintf(stderr,"\nUsage: makextal [-c A B C alpha beta gamma] [in.pdb \
[out.pdb]]\n");
   fprintf(stderr,"       -c Specify the unit cell parameters (default \
is to try to read\n");
   fprintf(stderr,"          them from the PDB file)\n");
   fprintf(stderr,"I/O is through standard input / standard output if no \
files are specified.\n");

   fprintf(stderr,"\nMake a crystal lattice from a PDB file. Places six \
molecules around the\n");
   fprintf(stderr,"central molecule which has the input \
coordinates.\n\n");
}


