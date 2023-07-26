/*************************************************************************

   Program:    
   File:       
   
   Version:    
   Date:       
   Function:   
   
   Copyright:  (c) UCL / Dr. Andrew C. R. Martin 2011
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

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bioplib/MathType.h"
#include "bioplib/pdb.h"
#include "bioplib/SysDefs.h"
#include "bioplib/macros.h"
#include "bioplib/angle.h"

/************************************************************************/
/* Defines and macros
*/

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
void Usage(void);


/************************************************************************/
int main(int argc, char **argv)
{
   int retval = 0;
   FILE *fp;
   
   
   if((argc==0) || !strncmp(argv[1], "-h", 2))
   {
      Usage();
   }
   else
   {
      if((fp=fopen(argv[1], "r"))==NULL)
      {
         fprintf(stderr, "Can't open file %s\n", argv[1]);
         retval = 1;
      }
      else
      {
         int natoms;
         PDB *pdb;
         if((pdb = ReadPDB(fp, &natoms))==NULL)
         {
            fprintf(stderr, "No atoms read from PDB file %s\n", argv[1]);
            retval = 1;
         }
         else
         {
            PDB *p;
            
            pdb = SelectCaPDB(pdb);
            for(p=pdb; p->next->next->next != NULL; NEXT(p))
            {
               REAL angle;
               PDB  *atom2 = p->next,
                    *atom3 = p->next->next,
                    *atom4 = p->next->next->next;
               char resid[16];
               
               angle = phi(p->x, p->y, p->z,
                           atom2->x, atom2->y, atom2->z, 
                           atom3->x, atom3->y, atom3->z, 
                           atom4->x, atom4->y, atom4->z);
               sprintf(resid, "%s%d%s", 
                       atom3->chain, atom3->resnum, atom3->insert);
               
               printf("%4s %-6s %8.3f\n",
                      atom3->resnam,
                      resid,
                      angle * 180.0 / PI);
            }
            
            FREELIST(pdb, PDB);
         }
      }
   }
   
   return(retval);
}

void Usage(void)
{
   printf("\npsudotorsion V0.1 (c) 2014 UCL, Dr. Andrew C.R. Martin   \n");

   printf("Usage: pseudotorsion pdbfile\n");

   printf("\nTakes a PDB file on the command line and calculates pseudo-torsion\n");
   printf("angles between the CA atoms. The torsions are assigned to the third\n");
   printf("CA atom in each group of four as done by Shirai et al.\n\n");
}
