/*************************************************************************

   Program:    pdbsphere
   File:       pdbsphere.c
   
   Version:    V1.1
   Date:       05.11.07
   Function:   Output all aminoacids within range from central aminoacid in
               a PDB file
   
   Copyright:  (c) UCL/Anja Baresic
   Author:     Anja Baresic
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
 
   Email:      anya@biochem.ucl.ac.uk
               
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
   Inputs ID of a central residue and a PDB file. 
   Outputs all residues from that PDB file, within 8 angstroms (default 
   range) of central residue's coordinates.

**************************************************************************

   Usage:
   ======
   
**************************************************************************

   Revision History:
   =================
   V1.0  23.10.07  Original
   V1.1  05.11.07  Added -r, -s, -h command line options
   V1.2  17.12.07  Changed output format to [chain]:resnum:[insert]
                   By: Anja

**************************************************************************/
/* Includes
*/
#include <stdio.h>
#include "bioplib/pdb.h"
#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/macros.h"

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
void FlagResiduesInRange(PDB *pdb, PDB *central, double radiusSq);
void WriteAtoms(PDB *pdb, FILE *out);
void WriteResidues(PDB *pdb, FILE *out);
BOOL ParseCmdLine(int argc, char **argv,char *resspec, char *InFile, 
                  char *OutFile, BOOL *summary, double *radiusSq);
void Usage(void);

/* ******************************************************************* */

/*>int main(int argc, char **argv)
   -------------------------------
   Input:  Central residue ID in [chain]num[insert] format . Input and 
           output filename can be specified after residue ID, if not, 
           program uses stdin and stdout.
   Output: Default:  Lists residues (in PDB format) within 8 angstroms 
                     range from central residue.
           -r radius:Changes range to radius
           -s       :Summary, outputs only a list of residues' IDs
           -h       :Prints out help
*/


int main(int argc, char **argv)
{
   FILE   *in  = stdin,
          *out = stdout;
   PDB    *pdb,
          *central;         
   int    natom;
   double radiusSq;
   char   resspec[MAXBUFF],
          InFile[MAXBUFF],
          OutFile[MAXBUFF];
    BOOL  summary;


   if (ParseCmdLine(argc, argv, resspec, InFile, OutFile, &summary, &radiusSq))
   {
      if (OpenStdFiles(InFile, OutFile, &in, &out))
      {
         if((pdb=ReadPDB(in, &natom))==NULL)
         {
            fprintf(stderr,"pdbsphere: No atoms read from PDB file\n");
            return(1);
         }
         
         else
         {        
            if ((central=FindResidueSpec(pdb, resspec))==NULL)
            {
               fprintf(stderr,"pdbsphere: Aminoacid %s not found in %s\n", resspec, InFile);
               return(1);
            }

            else                  
            {
              FlagResiduesInRange(pdb, central, radiusSq);

               if (summary)
               {
                  WriteResidues(pdb, out);
               }
               else
               {
                  WriteAtoms(pdb, out);
               }           
            }
         }
      }      
   }
   else 
   {
      Usage();
   }
   return 0;      
}





/**********************************************************************/
/* Subroutines
--------------
*/

/*> void FlagResiduesInRange(PDB *pdb, PDB *central, double *radiusSq)
-------------------------------------------------------------------
Input:   PDB    *pdb     
         PDB    *central    Pointer to the first atom of a central residue
         double radiusSq    To be flagged, atom has to be within that 
                            range (radius is squared for speed)
Output:  If any atom in a residue is within range from *central, marks all 
         atoms in that residue (sets occ parameter to 2.00) 
*/

void FlagResiduesInRange(PDB *pdb, PDB *central, double radiusSq)
{
   PDB *p,
       *q,
       *current,
       *nextPRes,
       *nextCurrentRes;
   BOOL aaInRange;
      

   nextPRes=FindNextResidue(central);
  
   for (current=pdb; current!=NULL; current=nextCurrentRes)
   {
     aaInRange=FALSE;
     nextCurrentRes=FindNextResidue(current);

     for (q=current; q!=nextCurrentRes; NEXT(q))
     {            
        for (p=central; p!=nextPRes; NEXT(p))
        {
           if (DISTSQ(p, q)<radiusSq)
           {
              aaInRange=TRUE;
              /*q=nextCurrentRes;
                break;*/
           }
        }
      }
     
        
        if (aaInRange)
        {
          for (q=current; q!=nextCurrentRes; NEXT(q))
          {
            q->occ=2.00;
          }    
        }
        
   }
} 




/************************************************************************/
/*>void WriteAtoms(PDB *pdb, FILE *out)
-----------------------------------------
Input:    PDB  *pdb   pointer to beginning of a linked list
          FILE *out   output file
Output    Writes atom's node in *out if it's occ is >1.90 (returns
          occ to original value)
*/

void WriteAtoms(PDB *pdb, FILE *out)
{  
   PDB *p;

   for (p=pdb; p!=NULL; NEXT(p))
   {
      if (p->occ > 1.90)
      {
         p->occ=1.00;
         WritePDBRecord(out, p);
      }
   }
}




/************************************************************************/
/*>void WriteResidues(PDB *pdb, FILE *out)
----------------------------------------------
Input:    PDB  *pdb   pointer to beginning of a linked list
          FILE *out   output file
Output    Writes a list of residues' IDs for residues marked with 
          occ>1.90 (returns occ to original value)
*/

void WriteResidues(PDB *pdb, FILE *out)
{

  PDB *p;
  
  for (p=pdb; p!=NULL; p=FindNextResidue(p))
   {
      if (p->occ > 1.90)
      {
         p->occ=1.00;
         fprintf(out, "%s:%d:%s\n", p->chain, p->resnum, p->insert);
      }
      
   }
  
}  



/***********************************************************************/
/*>void Usage(void)
--------------------
Writes help.
*/

void Usage(void)
{
   fprintf(stderr,"\n");
   fprintf(stderr,"PDBsphere V1.2 (c) Anja Baresic, UCL.\n");
   fprintf(stderr,"Last modified 17/12/07 by Anja Baresic, UCL.\n");
   fprintf(stderr,"\nUsage: \
PDBsphere [-s] [-r radius] [-h] resID [in.txt [out.txt]]\n");
   fprintf(stderr,"       -s  Output summary: only list of residue IDs.\n");
   fprintf(stderr,"       -r  Set your own allowed range to radius. \n");   
   fprintf(stderr,"Default behaviour is to output all atoms of \
a residue containing at least\n");
   fprintf(stderr,"one atom in range (from residue defined by resID), in a \
PDB format.\n");
   fprintf(stderr,"Default range is 8 angstroms.");
   fprintf(stderr,"\n");
   fprintf(stderr,"ResID is in form [c]num[i]where [c] is an optional chain \
specification, num\n");
   fprintf(stderr,"is a residue number and [i] is an optional insertion \
code.\n");
   fprintf(stderr,"\nPDBsphere writes all the residues within range of a \
residue with resID.\n");
   fprintf(stderr,"I/O through standard input/output if files not \
specified.\n\n");   
   
}



/**********************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *resspec, char *InFile, 
                    char *OutFile, BOOL *summary, double *radiusSq)
   ----------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *resspec     Central residue ID in [chain]num[insert]
                                format
            char   *InFile      Input file (or blank string)
            char   *OutFile     Output file (or blank string)
            BOOL   *summary     Should output be summarised?
                                (Default: no)
            double *radiusSq    Maximum allowed distance of residues' 
                                coordinates - squared
                                (Default:64, max range:8 angstroms)
   Returns: BOOL                Success?

   Parse the command line
   
   26.10.07 Original    By: Anya
*/
BOOL ParseCmdLine(int argc, char **argv,char *resspec, char *InFile, 
                  char *OutFile, BOOL *summary, double *radiusSq)
{
   argc--;
   argv++;

   InFile[0] = OutFile[0] = '\0';
   *summary = FALSE;
   *radiusSq = 64.00;
 
   if (argc<1)
   {      
     return (FALSE);
   }
      
   while(argc)
     {     
      if(argv[0][0] == '-')
      {
         if (argv [0][2]!='\0')
         {
           return(FALSE);
         }
         else
         {
            switch(argv[0][1])
            {
            case 'h':
               return(FALSE);
               break;
            case 's':
               *summary = TRUE;
               break;
            case 'r':
               argc--;
               argv++;
               if(argc < 0)
                  return(FALSE);
               if(!sscanf(argv[0],"%lf",radiusSq))
                  return(FALSE);
               else *radiusSq *= *radiusSq;  
               break;
            default:
               return(FALSE);
               break;
            }
         }
      }
      else
      {
         /* Check that there are 1, 2 or 3 arguments left               */
         if(argc<1 || argc > 3)
            return(FALSE);

         /* Copy the first to resspec                                    */
         if(argc)
         {
            strncpy(resspec, argv[0], MAXBUFF);
            argc--;
            argv++;
         }
         else
           return(FALSE);

         /* Copy the second to InFile                                    */
         if(argc)
         {
            strncpy(InFile, argv[0], MAXBUFF);
            argc--;
            argv++;
         }         

         /* If there's another, copy it to OutFile                       */
         if(argc)
         {
            strncpy(OutFile, argv[0], MAXBUFF);            
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
