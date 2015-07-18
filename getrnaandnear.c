/*************************************************************************

   Program:    
   File:       
   
   Version:    
   Date:       
   Function:   
   
   Copyright:  (c) UCL / Dr. Andrew C. R. Martin 2010
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   EMail:      martin@biochem.ucl.ac.uk
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
#include "bioplib/pdb.h"
#include "bioplib/seq.h"
#include "bioplib/macros.h"
#include "bioplib/SysDefs.h"

/************************************************************************/
/* Defines and macros
*/
/* Protein chain types                                                  */
#define TYPE_UNDEF   0
#define TYPE_RNA     1
#define TYPE_DNA     2
#define TYPE_PROTEIN 3
#define TYPE_HYBRID  4
#define TYPE_PEPTIDE 5
#define TYPE_NONSTD  6

/* Flags to indicate whether we will be keeping a chain                 */
#define CHAIN_NULL   0
#define CHAIN_RNA    1
#define CHAIN_NEAR   2

/* Distance cutoffs                                                     */
#define DISTCUTOFF 8
#define DISTCUTOFFSQ 64   /* 8A cutoff */

/* String buffer                                                        */
#define MAXBUFF 1024

/* Extra chain info                                                     */
typedef struct
{
   int type;
   REAL xmin, xmax, ymin, ymax, zmin, zmax;
} CHAININFO ;

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
int FindChainType(PDB *start, PDB *stop, BOOL doPep, BOOL doX);
BOOL CheckBounds(PDBCHAIN *chain1, PDBCHAIN *chain2);
void SetBounds(PDB *start, PDB *stop, REAL *xmin, REAL *xmax, REAL *ymin,
               REAL *ymax, REAL *zmin, REAL *zmax);
void AssignChainTypes(PDBSTRUCT *pdbs);
void FindChainsNearRNA(PDBSTRUCT *pdbs);
void PrintRNANearChains(FILE *out, PDBSTRUCT *pdbs);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile);
void Usage(void);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program

   03.03.10  Original   By: ACRM
*/
int main(int argc, char **argv)
{
   FILE *in = stdin;
   FILE *out = stdout;
   PDB  *pdb;
   PDBSTRUCT *pdbstruct;
   int  natoms;
   char infile[MAXBUFF],
        outfile[MAXBUFF];
   

   if(!ParseCmdLine(argc, argv, infile, outfile))
   {
      Usage();
      return(0);
   }
   
   /* Open files                                                        */
   if(!OpenStdFiles(infile, outfile, &in, &out))
   {
      fprintf(stderr, "Unable to open input or output files\n");
      return(1);
   }
   
   /* Read and format data                                              */
   if((pdb = ReadPDB(in, &natoms)) == NULL)
   {
      fprintf(stderr, "No atoms read from file %s\n", argv[0]);
      return(1);
   }
   if((pdbstruct = AllocPDBStructure(pdb)) == NULL)
   {
      fprintf(stderr, "No memory for PDB structure\n");
   }

   /* Process the data to find and print chains of interest             */
   AssignChainTypes(pdbstruct);
   FindChainsNearRNA(pdbstruct);
   PrintRNANearChains(out, pdbstruct);

   /* Clean up                                                          */
   FreePDBStructure(pdbstruct);
   FREELIST(pdb, PDB);

   return(0);
}


/************************************************************************/
/*>int FindChainType(PDB *start, PDB *stop, BOOL doPep, BOOL doX)
   --------------------------------------------------------------
   Input:   PDB  *start      Start of PDB linked list
            PDB  *stop       End of PDB linked list (or NULL)
            BOOL doPep       Count peptides as a separate type
            BOOL doX         Count chains with >25% of residues as non
                             standard amino acids as nonstandard chains
   Returns: int              Chain type

   Works out the type of a chain specified by the PDB pointer boundaries.
   The chain type is returned as 
   TYPE_UNDEF   0
   TYPE_RNA     1
   TYPE_DNA     2
   TYPE_PROTEIN 3
   TYPE_HYBRID  4
   TYPE_PEPTIDE 5 (if doPep true - otherwise TYPE_PROTEIN)
   TYPE_NONSTD  6 (if doX true - otherwise TYPE_PROTEIN)

   10.03.09 Original   By: ACRM
   30.03.09 Changed to use FindNextResidue() rather than checking every 
            atom. Also counts residues and checks for peptides and
            non-standard chains
*/
int FindChainType(PDB *start, PDB *stop, BOOL doPep, BOOL doX)
{
   PDB  *p;
   int  chaintype = TYPE_UNDEF;
   BOOL seenU = FALSE;
   int  rescount = 0;
   int  nonstd   = 0;
   
   for(p=start; (p!=stop && p!=NULL); p=FindNextResidue(p))
   {
      rescount++;
      if(throne(p->resnam) == 'X')
      {
         nonstd++;
      }
      
      /* If we don't have a chain type yet, and we see ACG then assume it
         is RNA
      */
      if(!strncmp(p->resnam, "  A", 3) ||
         !strncmp(p->resnam, "  C", 3) ||
         !strncmp(p->resnam, "  I", 3) ||
         !strncmp(p->resnam, "  G", 3))
      {
         if(chaintype != TYPE_RNA)
         {
            if(chaintype == TYPE_UNDEF)
            {
               chaintype = TYPE_RNA;
            }
            else if(chaintype == TYPE_PROTEIN)
            {
               chaintype = TYPE_HYBRID;
            }
         }
      }
      
      /* If we don't have a chain type yet, and we see U then it is RNA */
      else if(!strncmp(p->resnam, "  U", 3) ||
              !strncmp(p->resnam, " IU", 3))
      {
         seenU = TRUE;
         
         if(chaintype != TYPE_RNA)
         {
            if(chaintype == TYPE_UNDEF)
            {
               chaintype = TYPE_RNA;
            }
            else if((chaintype == TYPE_DNA) || 
                    (chaintype == TYPE_PROTEIN))
            {
               chaintype = TYPE_HYBRID;
            }
         }
      }
      
      /* If we see DA,DT,DC,DG or a T then it is DNA (even if previously 
         assigned as RNA)
      */
      else if(!strncmp(p->resnam, " DA", 3) ||
              !strncmp(p->resnam, " DT", 3) ||
              !strncmp(p->resnam, "  T", 3) ||
              !strncmp(p->resnam, " DC", 3) ||
              !strncmp(p->resnam, " DI", 3) ||
              !strncmp(p->resnam, " DG", 3))
      {
         if(chaintype != TYPE_DNA)
         {
            if((chaintype == TYPE_UNDEF) || (chaintype == TYPE_RNA))
            {
               if(seenU)
               {
                  chaintype = TYPE_HYBRID;
               }
               else
               {
                  chaintype = TYPE_DNA;
               }
            }
            else
            {
               chaintype = TYPE_HYBRID;
            }
         }
      }
      else   /* See anything else then it is protein                    */
      {
         if(chaintype != TYPE_PROTEIN)
         {
            if(chaintype == TYPE_UNDEF) 
            {
               chaintype = TYPE_PROTEIN;
            }
            else
            {
               chaintype = TYPE_HYBRID;
            }
         }
      }
   }

   if(doPep)
   {
      if((chaintype == TYPE_PROTEIN) && (rescount < 30))
         chaintype = TYPE_PEPTIDE;
   }
   if(doX)
   {
      if(((100*nonstd)/rescount) > 25)
         chaintype = TYPE_NONSTD;
   }

   return(chaintype);
}


/************************************************************************/
/*>void SetBounds(PDB *start, PDB *stop, REAL *xmin, REAL *xmax, 
                  REAL *ymin, REAL *ymax, REAL *zmin, REAL *zmax)
   --------------------------------------------------------------
   Run through a PDB linked list between start and stop and define the
   bounging box

   03.03.10  Original   By: ACRM
*/
void SetBounds(PDB *start, PDB *stop, REAL *xmin, REAL *xmax, REAL *ymin,
               REAL *ymax, REAL *zmin, REAL *zmax)
{
   PDB *p;

   *xmin = *xmax = start->x;
   *ymin = *ymax = start->y;
   *zmin = *zmax = start->z;

   for(p=start; p!=stop; NEXT(p))
   {
      if(p->x < *xmin)
      {
         *xmin = p->x;
      }
      else if(p->x > *xmax)
      {
         *xmax = p->x;
      }
      
      if(p->y < *ymin)
      {
         *ymin = p->y;
      }
      else if(p->y > *ymax)
      {
         *ymax = p->y;
      }
      
      if(p->z < *zmin)
      {
         *zmin = p->z;
      }
      else if(p->z > *zmax)
      {
         *zmax = p->z;
      }
   }
}

/************************************************************************/
/*>BOOL CheckBounds(PDBCHAIN *chain1, PDBCHAIN *chain2)
   ----------------------------------------------------
   Look to see if the bounding box of one chain is within DISTCUTOFF
   Angstroms of the bounding box of the second chain

   03.03.10  Original   By: ACRM
*/
BOOL CheckBounds(PDBCHAIN *chain1, PDBCHAIN *chain2)
{
   if(((((CHAININFO *)chain1->extras)->xmax) + DISTCUTOFF) >
      (((CHAININFO *)chain2->extras)->xmin))
      return(TRUE);
   if(((((CHAININFO *)chain1->extras)->ymax) + DISTCUTOFF) >
      (((CHAININFO *)chain2->extras)->ymin))
      return(TRUE);
   if(((((CHAININFO *)chain1->extras)->zmax) + DISTCUTOFF) >
      (((CHAININFO *)chain2->extras)->zmin))
      return(TRUE);
   
   if(((((CHAININFO *)chain2->extras)->xmax) + DISTCUTOFF) >
      (((CHAININFO *)chain1->extras)->xmin))
      return(TRUE);
   if(((((CHAININFO *)chain2->extras)->ymax) + DISTCUTOFF) >
      (((CHAININFO *)chain1->extras)->ymin))
      return(TRUE);
   if(((((CHAININFO *)chain2->extras)->zmax) + DISTCUTOFF) >
      (((CHAININFO *)chain1->extras)->zmin))
      return(TRUE);
   
   return(FALSE);
}


/************************************************************************/
/*>void AssignChainTypes(PDBSTRUCT *pdbs)
   --------------------------------------
   Runs through the chains, allocating memory for the 'extras' structure
   of type CHAININFO and setting the chain type, then calling the code
   to find the bounding box of each chain

   03.03.10  Original   By: ACRM
*/
void AssignChainTypes(PDBSTRUCT *pdbs)
{
   PDBCHAIN  *chain;
   int       chaintype;

   /* Assign a type to each chain                                       */
   for(chain=pdbs->chains; chain!=NULL; NEXT(chain))
   {
      chaintype = FindChainType(chain->start, 
                                ((chain->next)?chain->next->start:NULL),
                                FALSE, FALSE);
      if((chain->extras = (APTR *)malloc(sizeof(CHAININFO)))==NULL)
      {
         fprintf(stderr, "No memory for CHAININFO structure\n");
         exit(1);
      }

      if(chaintype == TYPE_RNA)
      {
         ((CHAININFO *)chain->extras)->type = CHAIN_RNA;
         fprintf(stderr,"RNA Chain %s\n", chain->chain);
      }
      else
      {
         ((CHAININFO *)chain->extras)->type = CHAIN_NULL;
      }

      /* Identify the chain boundaries                                  */
      SetBounds(chain->start, chain->stop, 
                &(((CHAININFO *)chain->extras)->xmin),
                &(((CHAININFO *)chain->extras)->xmax),
                &(((CHAININFO *)chain->extras)->ymin),
                &(((CHAININFO *)chain->extras)->ymax),
                &(((CHAININFO *)chain->extras)->zmin),
                &(((CHAININFO *)chain->extras)->zmax));
   }
}


/************************************************************************/
/*>void FindChainsNearRNA(PDBSTRUCT *pdbs)
   ---------------------------------------
   Does the main work of identifying chains which are near to RNA chains 

   03.03.10  Original   By: ACRM
*/
void FindChainsNearRNA(PDBSTRUCT *pdbs)
{
   PDBCHAIN  *chain1, *chain2;
   PDB       *p, *q;
      
   /* Run through each chain                                            */
   for(chain1=pdbs->chains; chain1!=NULL; NEXT(chain1))
   {
      /* If it's RNA                                                    */
      if(((CHAININFO *)chain1->extras)->type == CHAIN_RNA)
      {
         /* Run through all other chains                                */
         for(chain2=pdbs->chains; chain2!=NULL; NEXT(chain2))
         {
            /* If it's not RNA and not already flagged as near to RNA   */
            if((((CHAININFO *)chain2->extras)->type != CHAIN_RNA) &&
               (((CHAININFO *)chain2->extras)->type != CHAIN_NEAR))
            {
               /* If the bounding boxes of the chains are close enough  */
               if(CheckBounds(chain1, chain2))
               {
                  /* Run through atoms in first chain                   */
                  for(p=chain1->start; p!=chain1->stop; NEXT(p))
                  {
                     /* Run through atoms in second chain               */
                     for(q=chain2->start; q!=chain2->stop; NEXT(q))
                     {
                        /* If they are within specified distance        */
                        if(DISTSQ(p,q) < DISTCUTOFFSQ)
                        {
                           /* Flag the second chain as being near RNA   */
                           ((CHAININFO *)chain2->extras)->type = 
                              CHAIN_NEAR;
                           fprintf(stderr,"Near Chain %s\n", 
                                   chain2->chain);
                           /* Break out of the atom searches            */
                           goto breakout;
                        }
                     }  /* for() atoms in second chain                  */
                  }  /* for() atoms in first chain                      */
breakout:         continue;
               }  /* if boxes are close enough                          */
            }  /* if second chain not already RNA or NEAR               */
         }  /* for() each second chain                                  */
      }  /* if first chain is RNA                                       */
   }  /* for() each first chain                                         */
}


/************************************************************************/
/*>void PrintRNANearChains(FILE *out, PDBSTRUCT *pdbs)
   ---------------------------------------------------
   Print the RNA chains and the other chains which have been identified
   as being near to RNA chains

   03.03.10  Original   By: ACRM
*/
void PrintRNANearChains(FILE *out, PDBSTRUCT *pdbs)
{
   PDBCHAIN *chain;
   PDB      *p;

   /* Run through the chains                                            */
   for(chain=pdbs->chains; chain!=NULL; NEXT(chain))
   {
      /* If it's flagged as RNA or near to RNA, then print the records  */
      if((((CHAININFO *)chain->extras)->type == CHAIN_RNA) ||
         (((CHAININFO *)chain->extras)->type == CHAIN_NEAR))
      {
         for(p=chain->start; p!=chain->stop; NEXT(p))
         {
            WritePDBRecord(out, p);
         }
      }
   }
}

/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *ssfile, char *pdbfile, 
                     char *infile, char *outfile, char *chain, 
                     BOOL *verbose, int *flex)
   ----------------------------------------------------------------------
   Input:   int   argc     Command line argc
            char  **argv   Command line argv
   Output:  char  *infile  Input PDB filename (or blank string)
            char  *outfile Output PDB filename (or blank string)
   Returns: BOOL           Success

   Parse the command line

   03.03.10  Original   By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile)
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
         case 'h':
            return(FALSE);
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are 1 or 2 arguments left                  */
         if(argc > 2)
            return(FALSE);
         
         /* If another, Copy to infile                                  */
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
   Prints a usage message 

   03.03.10  Original   By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\ngetrnaandnear V1.0 (C) Dr. Andrew C.R Martin, \
UCL\n");

   fprintf(stderr,"\nUsage: getrnaandnear [infile [outfile]]\n");

   fprintf(stderr,"\ngetrnaandnear extracts RNA chains and chains that \
are within %.1f A\n", (REAL)DISTCUTOFF);
   fprintf(stderr,"of an RNA chain. It is used for cutting down huge \
structures such\n");
   fprintf(stderr,"as ribosomes before processing with programs like \
ligplot\n\n");
}
