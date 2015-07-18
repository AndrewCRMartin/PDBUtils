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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/array.h"
#include "bioplib/angle.h"



/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160

typedef struct _node
{
   struct _node *next,
                *prev;
   int          nextx,
                nexty,
                prevx,
                prevy,
                resnum;
   BOOL         done;
}  NODE;

typedef struct _cell
{
   int  count;
   NODE **nodes;
}  CELL;



/************************************************************************/
/* Globals
*/
CELL **gGrid = NULL;

/************************************************************************/
/* Prototypes
*/
int  main(int argc, char **argv);
BOOL AllocateGrid(int NGrid, int NStruc);
BOOL FillGrid(PDB *pdb, int StrucNum, int NGrid, 
              int *startx, int *starty);
BOOL DoKMRC(PDB **pdbs, int NStruc, int NGrid, int eta, 
            int *startx, int *starty);
int CalcBox(REAL angle, int NGrid);
void LinkPrevBox(int StrucNum, int xprev, int yprev, int xbox, int ybox);
BOOL SetThisBox(int StrucNum, int xbox, int ybox, int xprev, int yprev,
                int resnum);
void FindNextPoint(int i, int j, int resnum, int *inext, int *jnext);
BOOL CheckGridPoint(int i, int j, int NGrid, int NStruc, int eta);
BOOL Consecutive(int i, int j, int resnum, int NStruc, int i1, int j1, 
                 int *group1, int *group2);
   

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   03.09.97 Original   By: ACRM
*/
int main(int argc, char **argv)
{
   int  NGrid,
        eta,
        i,
        natoms,
        NStruc,
        *startx,
        *starty;
   PDB  **pdbs;
   char buffer[MAXBUFF];
   FILE *fp;
   
   
   /* Read the number of divisions along phi and psi                    */
   PROMPT(stdin, "Enter number of grid points: ");
   fgets(buffer,MAXBUFF,stdin);
   TERMINATE(buffer);
   sscanf(buffer,"%d", &NGrid);
   
   /* Read the tolerence (number of grid blocks grouped as 1)           */
   PROMPT(stdin, "Enter tolerence (in grid units): ");
   fgets(buffer,MAXBUFF,stdin);
   TERMINATE(buffer);
   sscanf(buffer,"%d", &eta);
   
   /* Read number of structures to fit                                  */
   PROMPT(stdin, "Enter number of structures to compare: ");
   fgets(buffer,MAXBUFF,stdin);
   TERMINATE(buffer);
   sscanf(buffer,"%d", &NStruc);

   /* Allocate memory for the grid                                      */
   if(!AllocateGrid(NGrid, NStruc))
   {
      fprintf(stderr,"Unable to allocate grid\n");
      return(1);
   }

   /* Allocate memory for the PDB linked lists and starts               */
   if((pdbs = (PDB **)malloc(NStruc * sizeof(PDB *)))==NULL)
   {
      fprintf(stderr,"No memory for array of PDB pointers\n");
      return(1);
   }
   if((startx = (int *)malloc(NStruc * sizeof(int)))==NULL)
   {
      fprintf(stderr,"No memory for startx array\n");
      return(1);
   }
   if((starty = (int *)malloc(NStruc * sizeof(int)))==NULL)
   {
      fprintf(stderr,"No memory for starty array\n");
      return(1);
   }
   
   
   /* Read the structures and fill in the grid                          */
   for(i=0; i<NStruc; i++)
   {
      PROMPT(stdin, "Enter structure name: ");
      fgets(buffer,MAXBUFF,stdin);
      TERMINATE(buffer);

      if((fp=fopen(buffer,"r"))==NULL)
      {
         fprintf(stderr,"Unable to read PDB file: %s\n",buffer);
         return(1);
      }
      if((pdbs[i] = ReadPDB(fp, &natoms))==NULL)
      {
         fprintf(stderr,"No atoms read from PDB file: %s\n",buffer);
         return(1);
      }
      
      FillGrid(pdbs[i], i, NGrid, &(startx[i]), &(starty[i]));
      fclose(fp);
   }

   /* Now run the actual KMRC algorithm                                 */
   if(!DoKMRC(pdbs, NStruc, NGrid, eta, startx, starty))
   {
      fprintf(stderr,"No memory for KMRC algorithm\n");
      return(1);
   }

   return(0);
}


/************************************************************************/
/*>BOOL AllocateGrid(int NGrid, int NStruc)
   ----------------------------------------
   03.09.97 Original   By: ACRM
*/
BOOL AllocateGrid(int NGrid, int NStruc)
{
   int i, j, k;
   
   /* Create the 2D array                                               */
   if((gGrid = (CELL **)Array2D(sizeof(CELL), NGrid, NGrid))==NULL)
      return(FALSE);

   /* Now 0 the counts in the grid and allocate arrays and NULL them    */
   for(i=0; i<NGrid; i++)
   {
      for(j=0; j<NGrid; j++)
      {
         gGrid[i][j].count = 0;
         if((gGrid[i][j].nodes = (NODE **)malloc(NStruc * sizeof(NODE *)))
            == NULL)
            return(FALSE);
         for(k=0; k<NStruc; k++)
            gGrid[i][j].nodes[k] = NULL;
      }
   }
   
   return(TRUE);
}

/************************************************************************/
/*>BOOL FillGrid(PDB *fullpdb, int StrucNum, int NGrid, 
                 int *startx, int *starty)
   -----------------------------------------------------
   04.09.97 Original   By: ACRM
*/
BOOL FillGrid(PDB *fullpdb, int StrucNum, int NGrid, 
              int *startx, int *starty)
{
   PDB  *p, *pdb,
        *p1, *p2, *p3, *p4;
   char *sel[4];
   int  natoms, resnum,
        xbox,  ybox,
        xprev, yprev;
   REAL Phi, Psi, Omega;
   
   
   SELECT(sel[0],"CA  ");
   SELECT(sel[1],"N   ");
   SELECT(sel[2],"C   ");

   xprev = yprev = (-1);
   xbox  = ybox  = (-1); /* Indicate start of structure                 */

   if((pdb = SelectAtomsPDB(fullpdb,3,sel,&natoms))==NULL)
   {
      fprintf(stderr,"Unable to select backbone atoms from PDB \
file (no memory?)\n");
      return(FALSE);
   }
   
   /* Walk the linked list and calculate torsions                       */
   Phi     = Psi     = Omega   = 9999.0;
   p1      = p2      = p3      = p4     = NULL;
   resnum  = 2;
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(!strncmp(p->atnam,"C   ",4) &&
         (Phi != 9999.0)             &&
         (Psi != 9999.0))
      {
         xprev = xbox;
         yprev = ybox;
         xbox = CalcBox(Phi, NGrid);
         ybox = CalcBox(Psi, NGrid);
         if(xprev == (-1))
         {
            *startx = xbox;
            *starty = ybox;
         }
         LinkPrevBox(StrucNum, xprev, yprev, xbox, ybox);
         if(!SetThisBox(StrucNum, xbox, ybox, xprev, yprev, resnum++))
            return(FALSE);
      }
         
      /* Get pointers to four atoms in sequence                         */
      p1 = p;
      p2 = p->next;
      if(p2 != NULL) p3 = p2->next;
      if(p3 != NULL) p4 = p3->next;

      if(p2==NULL || p3==NULL || p4==NULL)
         break;
         
      if(!strncmp(p->atnam,"N   ",4))
         Psi = phi(p1->x, p1->y, p1->z,
                   p2->x, p2->y, p2->z,
                   p3->x, p3->y, p3->z,
                   p4->x, p4->y, p4->z);
      else if(!strncmp(p->atnam,"CA  ",4))
         Omega = phi(p1->x, p1->y, p1->z,
                     p2->x, p2->y, p2->z,
                     p3->x, p3->y, p3->z,
                     p4->x, p4->y, p4->z);
      else if(!strncmp(p->atnam,"C   ",4))
         Phi = phi(p1->x, p1->y, p1->z,
                   p2->x, p2->y, p2->z,
                   p3->x, p3->y, p3->z,
                   p4->x, p4->y, p4->z);
   }
   
   return(TRUE);
}

/************************************************************************/
/*>int CalcBox(REAL angle, int NGrid)
   ----------------------------------
   04.09.97 Original   By: ACRM
*/
int CalcBox(REAL angle, int NGrid)
{
   return((int)((simpleangle(angle) * (REAL)NGrid)/((REAL)2.0 * PI)));
}

/************************************************************************/
/*>void LinkPrevBox(int StrucNum, int xprev, int yprev, int xbox, 
                    int ybox)
   --------------------------------------------------------------
   04.09.97 Original   By: ACRM
*/
void LinkPrevBox(int StrucNum, int xprev, int yprev, int xbox, int ybox)
{
   CELL *cell;
   NODE *node;
   
   if((xprev != (-1)) && (yprev != (-1)))
   {
      /* Just set these pointers to make things easier                  */
      cell = &(gGrid[xprev][yprev]);
      node = cell->nodes[StrucNum];
      LAST(node);
      
      node->nextx = xbox;
      node->nexty = ybox;
   }
}

/************************************************************************/
/*>BOOL SetThisBox(int StrucNum, int xbox, int ybox, int xprev, int yprev,
                   int resnum)
   -----------------------------------------------------------------------
   04.09.97 Original   By: ACRM
*/
BOOL SetThisBox(int StrucNum, int xbox, int ybox, int xprev, int yprev,
                int resnum)
{
   CELL *cell;
   NODE *node;
   
   /* Just set these pointers to make things easier                     */
   cell = &(gGrid[xbox][ybox]);
   node = cell->nodes[StrucNum];

   /* Allocate space for another node                                   */
   if(node==NULL)
   {
      INITPREV(cell->nodes[StrucNum], NODE);
      node = cell->nodes[StrucNum];
   }
   else
   {
      LAST(node);
      ALLOCNEXTPREV(node, NODE);
   }
   if(node == NULL)
      return(FALSE);

   /* Fill in this node                                                 */
   node->nextx  = (-1);
   node->nexty  = (-1);
   node->prevx  = xprev;
   node->prevy  = yprev;
   node->done   = FALSE;
   node->resnum = resnum;

   return(TRUE);
}

/************************************************************************/
/*>void FindNextPoint(int i, int j, int resnum, int *inext, int *jnext)
   --------------------------------------------------------------------
   For structure 1 finds the next point from i,j for resnum

   04.09.97 Original   By: ACRM
*/
void FindNextPoint(int i, int j, int resnum, int *inext, int *jnext)
{
   CELL *cell;
   NODE *node;

   /* Just set these pointers to make things easier                     */
   cell = &(gGrid[i][j]);
   for(node = cell->nodes[0]; node!=NULL; NEXT(node))
   {
      if(node->resnum == resnum)
      {
         *inext = node->nextx;
         *jnext = node->nexty;
         return;
      }
   }

   *inext = (-1);
   *jnext = (-1);
}

/************************************************************************/
/*>BOOL CheckGridPoint(int x, int y, int NGrid, int NStruc, int eta)
   ------------------------------------------------------------------
   Tests whether this grid point has hits in all NStruc structures

   04.09.97 Original   By: ACRM
*/
BOOL CheckGridPoint(int x, int y, int NGrid, int NStruc, int eta)
{
   CELL *cell;
   int  i,j,k,i1,j1;
   BOOL *ok, retval;

   ok = (BOOL *)malloc(NStruc*sizeof(BOOL));

   ok[0] = TRUE;
   /* Run through the structures, assuming that they are not found in
      this cell group
   */
   for(k=1; k<NStruc; k++)
   {
      ok[k] = FALSE;

      /* Run through the cells in this group setting i1 and j1 to account
         for wrap-around from one side of the grid to the other
      */
      for(i=(x-eta); i<=(x+eta); i++)
      {
         if(i<0) 
            i1=NGrid+i;
         else if(i>=NGrid) 
            i1=i-NGrid;
         else
            i1=1;
         
         for(j=(y-eta); j<=(y+eta); j++)
         {
            if(j<0) 
               j1=NGrid+j;
            else if(j>=NGrid) 
               j1=j-NGrid;
            else
               j1=1;

            /* If this cell has an entry for structure k, then set the
               flag to TRUE to say this structure is represented in this
               cell group
            */
            cell = &(gGrid[i1][j1]);
            if(cell->nodes[k] != NULL)
               ok[k] = TRUE;
         }
      }
   }

   /* Now run through the set of flags and see whether all structures were
      represented in this group
   */
   retval = TRUE;
   for(k=0; k<NStruc; k++)
   {
      if(!ok[k])
         retval = FALSE;
   }
   free(ok);

   return(retval);
}

/************************************************************************/
/*>BOOL DoKMRC(PDB **pdbs, int NStruc, int NGrid, int eta, 
               int *startx, int *starty)
   --------------------------------------------------------
*/
BOOL DoKMRC(PDB **pdbs, int NStruc, int NGrid, int eta, 
            int *startx, int *starty)
{
   int  i,  j, 
        resnum,
        inext, 
        jnext,
        *group1,
        *group2;

   /* Allocate memory for the next and prev arrays                      */
   if((group1 = (int *)malloc(NStruc * sizeof(int)))==NULL)
      return(FALSE);
   if((group2 = (int *)malloc(NStruc * sizeof(int)))==NULL)
      return(FALSE);

   /* Run along structure 1                                             */
   i = startx[0];
   j = starty[0];
   resnum = 2;     /* Our first residue is number 2                     */
   
   /* If this point has equivalents in the other structures             */
   if(CheckGridPoint(i, j, NGrid, NStruc, eta))
   {
      FindNextPoint(i, j, resnum, &inext, &jnext);
      if(CheckGridPoint(inext, jnext, NGrid, NStruc, eta))
      {
         if(Consecutive(i,j, resnum, NStruc, inext,jnext, group1,group2))
         {
            /* We have a 2 residue fragment which is common, so we
               stash these residue numbers as part of the common set
HERE
            */


            /* Step along the NStruc structures a residue at a time
               adding them to the common set if they are in equivalent
               positions
HERE
            */


         }
      }
   }

   free(group1);
   free(group2);
   return(TRUE);
}

/************************************************************************/
/*>BOOL Consecutive(int i1, int j1, int resnum, int NStruc, 
                    int i2, int j2, int *group1, int *group2)
   ----------------------------------------------------------
   Given that point i1,j1 and point i2,j2 both have all structures
   represented, tests whether they have consecutive residue numbers in
   each structure.
   If so, fills in group1 and group2 with the residue numbers from the
   structures at the two grid points.

TODO: Somehow this needs to handle structures which have multiple
   2-residue fragments where we will return multiple groups (i.e.
   group{1|2} need to be linked lists of malloc'd groups.
*/
BOOL Consecutive(int i1, int j1, int resnum, int NStruc, int i2, int j2, 
                 int *group1, int *group2)
{
   int k;
   CELL *cell1,
        *cell2;
   NODE *node1;

   /* We know this one...                                               */
   group1[0] = resnum;

   /* Just set these pointers to make things easier                     */
   cell1 = &(gGrid[i1][j1]);
   cell2 = &(gGrid[i2][j2]);

   for(k=1; k<NStruc; k++)
   {
      for(node1 = cell1->nodes[0]; node1!=NULL; NEXT(node1))
      {
         if(node1->resnum == resnum)
         {
            break;
         }
      }
      if(node1==NULL)
      {
         fprintf(stderr,"Consecutive(): Internal confusion!\n");
         exit(1);
      }
   }
   
   
   return(TRUE);
}

