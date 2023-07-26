/*************************************************************************

   Program:    
   File:       
   
   Version:    
   Date:       
   Function:   
   
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

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/pdb.h"

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

/************************************************************************/
int main(int argc, char **argv)
{
   FILE *fp, *pdbfp;
   char buffer[MAXBUFF],
        filename[MAXBUFF],
        type[MAXBUFF];
   PDB *pdb;
   REAL resol, RFac;
   int natoms, StrucType,
       ModelCount = 0,
       XtalCount  = 0,
       NMRCount   = 0,
       UnknownCount = 0;
   
   if((pdbfp = fopen(argv[1],"r"))!=NULL)
   {
      if(GetResolPDB(pdbfp,&resol,&RFac,&StrucType))
      {
         switch(StrucType)
         {
         case STRUCTURE_TYPE_XTAL:
            strcpy(type,"crystal");
            XtalCount++;
            break;
         case STRUCTURE_TYPE_NMR:
            strcpy(type,"NMR");
            NMRCount++;
            break;
         case STRUCTURE_TYPE_MODEL:
            strcpy(type,"model");
            ModelCount++;
            break;
         default:
            strcpy(type,"unknown");
            UnknownCount++;
            break;
         }
         printf("%.2f\n",resol);
      }
      else
      {
         printf("%s: No valid info\n",buffer);
      }
      
      fclose(pdbfp);
   }
   else
   {
      fprintf(stderr,"Unable to open file %s\n",argv[1]);
      return(1);
   }
   return(0);
}

               
         
      
      
  
