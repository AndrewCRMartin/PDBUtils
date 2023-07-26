/*************************************************************************

   Program:    splitpdb
   File:       splitpdb.c
   
   Version:    V1.0
   Date:       22.10.96
   Function:   Splits a PDB file using REMARK fields
   
   Copyright:  (c) Dr. Andrew C. R. Martin 1996
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   Phone:      (Home) +44 (0)1372 275775
               (Work) +44 (0)171 419 3890
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
#include <string.h>

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
   FILE *fp,
        *out = NULL;
   char buffer[MAXBUFF],
        buff2[MAXBUFF],
        *cnam,
        *chp;
   
   if(argc != 2)
   {
      fprintf(stderr,"Usage: split input.pdb\n");
      return(0);
   }
   
   if((fp=fopen(argv[1],"r"))==NULL)
   {
      fprintf(stderr,"Unable to open file: %s\n",argv[1]);
      return(1);
   }
   
   while(fgets(buffer, MAXBUFF, fp))
   {
      if(!strncmp(buffer,"REMARK",6))
      {
         /* Make a copy                                                 */
         strcpy(buff2, buffer);
         
         /* Get the codename out of this record                         */
         cnam = buffer+7;
         while(*cnam==' ') cnam++;
         for(chp=cnam; *chp && *chp!=' '; chp++);
         *chp='\0';

         /* Close any already-open file                                 */
         if(out != NULL)
         {
            fclose(out);
            out = NULL;
         }
         
         /* Open the new file                                           */
         if((out=fopen(cnam,"w"))==NULL)
         {
            fprintf(stderr,"Unable to open file for writing: %s\n",cnam);
         }
         else
         {
            fputs(buff2,out);
         }
      }
      else if(!strncmp(buffer, "ATOM  ",6))
      {
         if(out!=NULL)
         {
            fputs(buffer,out);
         }
      }
   }
   if(out!=NULL)
   {
      fclose(out);
   }
   return(0);
}

            
         

         
