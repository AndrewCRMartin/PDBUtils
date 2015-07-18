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
   EMail:      andrew@bioinf.org.uk
               andrew.martin@ucl.ac.uk
               
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
   FILE *pdbfp;
   char buffer[MAXBUFF],
        type[MAXBUFF];
   REAL resol, RFac;
   int StrucType;
   
   if((pdbfp = fopen(argv[1],"r"))!=NULL)
   {
      if(blGetResolPDB(pdbfp,&resol,&RFac,&StrucType))
      {
         switch(StrucType)
         {
         case STRUCTURE_TYPE_XTAL:
            strcpy(type,"crystal");
            break;
         case STRUCTURE_TYPE_NMR:
            strcpy(type,"NMR");
            break;
         case STRUCTURE_TYPE_MODEL:
            strcpy(type,"model");
            break;
         case STRUCTURE_TYPE_ELECTDIFF:
            strcpy(type,"ElectronDiffraction");
            break;
         case STRUCTURE_TYPE_FIBER:
            strcpy(type,"FiberDiffraction");
            break;
         case STRUCTURE_TYPE_SSNMR:
            strcpy(type,"SolidStateNMR");
            break;
         case STRUCTURE_TYPE_NEUTRON:
            strcpy(type,"NeutronScattering");
            break;
         case STRUCTURE_TYPE_EM:
            strcpy(type,"ElectronMicroscopy");
            break;
         case STRUCTURE_TYPE_SOLSCAT:
            strcpy(type,"SolutionScattering");
            break;
         case STRUCTURE_TYPE_IR:
            strcpy(type,"InfraredSpectroscopy");
            break;
         case STRUCTURE_TYPE_POWDER:
            strcpy(type,"PowderDiffraction");
            break;
         case STRUCTURE_TYPE_FRET:
            strcpy(type,"FlourescenceTransfer");
            break;
         default:
            strcpy(type,"unknown");
            break;
         }
         printf("%s, %.2fA/%.2f%%\n",type,resol,RFac*100.0);
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

               
         
      
      
  
