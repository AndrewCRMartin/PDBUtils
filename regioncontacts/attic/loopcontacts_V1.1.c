/*************************************************************************

   Program:    loopcontacts
   File:       loopcontacts.c
   
   Version:    V1.1
   Date:       01.10.96
   Function:   Analyse contacts between regions of structure (typically
               loops).
   
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
   Calculates loop takeoff angles.

   The algorithm we use is to calculate the CofG of the loop, the midpoint
   of the N-C-terminal vector and the CofG of the (5) residues on either
   side of the loop. The midpoint is moved to the origin and the loop
   is rotated such that the C-terminus is along the +ve x-axis and the
   framework CofG is on the xy-plane with -ve y.

   Two angles are then provided; one is the in-plane angle which describes
   how skewed the loop is towards the N- or C-terminus. This angle is
   +ve for skew towards the C-terminus.

   The other is the out-of-plane angle which describes how the loop
   flaps up and down with respect to the plane of the adjoining secondary
   structure. With the N-ter to your left and the C-ter to your right
   with the loop at the top, if it flops towards you, the angle is +ve.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  30.09.96 Original
   V1.1  01.10.96 Added analysis of contacting residues rather than
                  contacting pairs

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bioplib/MathType.h"
#include "bioplib/MathUtil.h"
#include "bioplib/SysDefs.h"
#include "bioplib/pdb.h"
#include "bioplib/general.h"
#include "bioplib/parse.h"
#include "bioplib/macros.h"
#include "bioplib/matrix.h"
#include "bioplib/angle.h"


/************************************************************************/
/* Defines and macros
*/
#define KEY_LOOP             0
#define KEY_QUIT             1
#define KEY_PDB              2
#define KEY_CONTACTS         3
#define KEY_RESIDUES         4
#define KEY_EXIT             5
#define PARSER_NCOMM         6
#define PARSER_MAXSTRPARAM   3
#define PARSER_MAXSTRLEN     80
#define PARSER_MAXREALPARAM  1

#define MAXBUFF              160
#define DISTMAX              ((REAL)4.0)

typedef struct _region
{
   struct _region *next;
   PDB  *start,
        *end;
   char startres[16],
        endres[16];
}  REGION;

typedef struct _contact
{
   struct _contact *next;
   char res1[16],
        res2[16];
   BOOL done;
   int  count;
}  CONTACT;

typedef struct _reslist
{
   struct _reslist *next;
   char res[16];
   int count;
}  RESLIST;


/************************************************************************/
/* Globals
*/
static MKeyWd sKeyWords[PARSER_NCOMM];         /* Parser keywords       */
static char   *sStrParam[PARSER_MAXSTRPARAM];  /* Parser string params  */
static REAL   sRealParam[PARSER_MAXREALPARAM], /* Parser real params    */
              gDistSqMax = DISTMAX * DISTMAX;
static int     gLoopCount = 0;
static REGION  *gLoops    = NULL;
static CONTACT *gContacts = NULL,
               *gTotal    = NULL;
static RESLIST *gResList  = NULL,
               *gCResList = NULL,
               *gCRTotal  = NULL;


/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile);
BOOL SetupParser(void);
BOOL ProcessInFile(FILE *in, FILE *out);
BOOL StoreLoop(char *startres, char *endres);
void PatchRegions(PDB *pdb);
BOOL ProcessFile(FILE *out, char *filename);
BOOL UpdateResCounts(void);
BOOL DoContactAnalysis(FILE *out);
void ClearContacts(void);
BOOL StoreContact(PDB *p, PDB *q);
void DisplayContacts(FILE *out, CONTACT *clist, REAL cutoff);
void DisplayReslist(FILE *out, RESLIST *rlist, REAL cutoff);
BOOL UpdateTotals(void);
CONTACT *GotContact(char *resspec1, char *resspec2, CONTACT *clist);
void Usage(void);
RESLIST *InResList(char *rspec, RESLIST *rlist);
int CountRes(char *rspec);
char *PDBResSpec(PDB *p);
BOOL StoreCRes(PDB *p);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for counting interloop contacts

   25.09.96 Original   By: ACRM
*/
int main(int argc, char **argv)
{
   char infile[MAXBUFF],
        outfile[MAXBUFF];
   FILE *in = stdin,
        *out = stdout;


   if(ParseCmdLine(argc, argv, infile, outfile))
   {
      if(SetupParser())
      {
         if(OpenStdFiles(infile, outfile, &in, &out))
         {
            if(ProcessInFile(in, out))
               return(0);
         }
      }
      return(1);
   }
   else
   {
      Usage();
   }
   return(0);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
   Returns: BOOL                Success?

   Parse the command line
   
   25.09.96 Original    By: ACRM
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
/*>BOOL SetupParser(void)
   ----------------------
   Returns:   BOOL         Success of memory allocations

   Set up the command parser

   25.09.96 Original    By: ACRM
*/
BOOL SetupParser(void)
{
   int i;

   
   /* Allocate memory for the string parameters                         */
   for(i=0; i<PARSER_MAXSTRPARAM; i++)
   {
      if((sStrParam[i] = 
         (char *)malloc(PARSER_MAXSTRLEN * sizeof(char)))==NULL)
      {
         int j;
         
         for(j=0;j<i;j++) 
            free(sStrParam[j]);
         fprintf(stderr,"No memory for parser string array\n");
         
         return(FALSE);
      }
   }
   
   /* Set up the keywords                                               */
   MAKEMKEY(sKeyWords[KEY_LOOP],       "LOOP",       STRING, 2, 2);
   MAKEMKEY(sKeyWords[KEY_EXIT],       "EXIT",       NUMBER, 0, 0);
   MAKEMKEY(sKeyWords[KEY_QUIT],       "QUIT",       NUMBER, 0, 0);
   MAKEMKEY(sKeyWords[KEY_PDB],        "PDB",        STRING, 1, 1);
   MAKEMKEY(sKeyWords[KEY_CONTACTS],   "CONTACTS",   NUMBER, 1, 1);
   MAKEMKEY(sKeyWords[KEY_RESIDUES],   "RESIDUES",   NUMBER, 1, 1);
   
   /* Check all allocations OK                                          */
   for(i=0; i<PARSER_NCOMM; i++)
   {
      if(sKeyWords[i].name == NULL)
      {
         fprintf(stderr,"No memory for keywords, or keyword undefined\n");
         return(FALSE);
      }
   }
   
   return(TRUE);
}


/************************************************************************/
/*>BOOL ProcessInFile(FILE *in, FILE *out)
   ---------------------------------------
   Input:   FILE  *in          Input file pointer
            FILE  *out         Output file pointer
   Returns: BOOL               Success? Fails if illegal input encountered

   Main loop to handle the command parser for the control file

   25.09.96 Original   By: ACRM
   01.10.96 Removed unused variable
*/
BOOL ProcessInFile(FILE *in, FILE *out)
{
   char buffer[MAXBUFF];
   int  NParams,
        key;


   PROMPT(in,"LoopContacts> ");
   
   while(fgets(buffer,MAXBUFF,in))
   {
      TERMINATE(buffer);

      if(!upstrncmp(buffer,"ECHO",4))
      {
         fprintf(out,"%s\n",buffer+5);
      }
      else
      {
         key = mparse(buffer, PARSER_NCOMM, sKeyWords, sRealParam,
                      sStrParam, &NParams);

         switch(key)
         {
         case PARSE_COMMENT:
            break;
         case PARSE_ERRC:
            fprintf(stderr,"Error in command: %s\n",buffer);
            break;
         case PARSE_ERRP:
            fprintf(stderr,"Error in parameters: %s\n",buffer);
            break;
         case KEY_LOOP:
            if(!StoreLoop(sStrParam[0], sStrParam[1]))
            {
               fprintf(stderr,"No memory for loop data!\n");
               return(FALSE);
            }
            break;
         case KEY_CONTACTS:
            fprintf(out,"----------------------------------------------\
----------------------------\n");
            fprintf(out,"Contacts made in >%.2f of cases\n",
                    sRealParam[0]);
            fprintf(out,"===============================\n");
            
            DisplayContacts(out, gTotal, sRealParam[0]);
            break;
         case KEY_RESIDUES:
            fprintf(out,"----------------------------------------------\
----------------------------\n");
            fprintf(out,"Residues making contacts in >%.2f of cases\n",
                    sRealParam[0]);
            fprintf(out,"==========================================\n");
            
            DisplayReslist(out, gCRTotal, sRealParam[0]);
            break;
         case KEY_EXIT:
            fprintf(out,"\n--------------------------------------------\
------------------------------\n");
            fprintf(out,"Summary Contacts\n================\n");
            DisplayContacts(out, gTotal, (REAL)(-1.0));
            fprintf(out,"\n--------------------------------------------\
------------------------------\n");
            fprintf(out,"Summary Contacting Residues\n=================\
==========\n");
            DisplayReslist(out, gCRTotal, (REAL)(-1.0));
            return(TRUE);
            break;
         case KEY_QUIT:
            return(TRUE);
            break;
         case KEY_PDB:
            if(!ProcessFile(out, sStrParam[0]))
            {
               fprintf(stderr,"Error processing PDB file: %s\n",
                       sStrParam[0]);
            }
            break;
         default:
            break;
         }
      }
      PROMPT(in,"LoopContacts> ");
   }  

   return(TRUE);
}


/************************************************************************/
/*>BOOL StoreLoop(char *startres, char *endres)
   --------------------------------------------
   Stores details of a loop specified in the control file

   25.09.96 Original   By: ACRM
*/
BOOL StoreLoop(char *startres, char *endres)
{
   static REGION *p;

   
   if(gLoops == NULL)
   {
      INIT(gLoops, REGION);
      p = gLoops;
   }
   else
   {
      ALLOCNEXT(p, REGION);
   }

   if(p==NULL)
   {
      FREELIST(gLoops, REGION);
      return(FALSE);
   }
   
   strncpy(p->startres, startres, 16);
   strncpy(p->endres,   endres,   16);
   
   return(TRUE);
}


/************************************************************************/
/*>void PatchRegions(PDB *pdb)
   ---------------------------
   Patches the actual PDB pointers into the linked list of loop regions.

   25.09.96 Original   By: ACRM
*/
void PatchRegions(PDB *pdb)
{
   REGION *r;
   PDB    *p;
   BOOL   found;

   
   /* Go through each of the loop specifications                        */
   for(r=gLoops; r!=NULL; NEXT(r))
   {
      /* Set the PDB pointers to NULLs                                  */
      r->start = r->end = NULL;
      found = FALSE;

      /* Run through the PDB linked list                                */
      for(p=pdb; p!=NULL; NEXT(p))
      {
         /* If we're in range                                           */
         if(InPDBZoneSpec(p, r->startres, r->endres))
         {
            /* Store start of zone                                      */
            if(r->start == NULL)
               r->start = p;

            /* Update end of zone                                       */
            r->end = p;
            found  = TRUE;
         }
         else
         {
            /* For efficiency, if we've found anything in zone before then
               we've come to the end of the zone and won't find anything
               else
            */
            if(found)
               break;
         }
      }  /* Loop through PDB linked list                                */

      /* Step the end of zone on so that we can use a simple not-equal
         test
      */
      if(r->end != NULL)
         r->end = r->end->next;
      
   }  /* Loop through ranges                                            */
}


/************************************************************************/
/*>BOOL ProcessFile(FILE *out, char *filename)
   -------------------------------------------
   Process a PDB file. The loop zone specifications must already have
   been given. Opens and reads the specified PDB file. Calls code to
   patch the PDB pointers into the linked list of loop specifications.
   Calls code to update the counts of each residue number found within
   the loop regions then does the contact analysis and frees the PDB
   linked list.

   25.09.96 Original   By: ACRM
*/
BOOL ProcessFile(FILE *out, char *filename)
{
   FILE *fp;
   PDB  *pdb;
   BOOL retval = TRUE;
   int  natoms;

   
   /* Open and read PDB file                                            */
   if((fp=fopen(filename,"r"))==NULL)
   {
      fprintf(stderr,"Unable to open file\n");
      return(FALSE);
   }
   if((pdb=ReadPDB(fp, &natoms))==NULL)
   {
      fprintf(stderr,"No atoms read\n");
      fclose(fp);
      return(FALSE);
   }
   fclose(fp);

   gLoopCount++;

   /* Patch the PDB pointers into the linked list of loops              */
   PatchRegions(pdb);
   
   /* Update the counts for each residue number so we can calculate
      percentages later
   */
   UpdateResCounts();

   /* Perform the contact analysis                                      */
   if(!DoContactAnalysis(out))
      retval = FALSE;

   /* Free the PDB linked list                                          */
   FREELIST(pdb, PDB);
   
   return(retval);
}


/************************************************************************/
/*>BOOL UpdateResCounts(void)
   --------------------------
   Maintains a linked list of counts of each residue number observed in
   the loop zones. This is used to calculate percentages for those 
   residues involved in contacts.

   30.09.96 Original   By: ACRM
   01.10.96 Changed call to InResList()
*/
BOOL UpdateResCounts(void)
{
   static RESLIST *sl;
   RESLIST        *l;
   REGION         *r;
   PDB            *p;
   
   
   /* For each region                                                   */
   for(r=gLoops; r!=NULL; NEXT(r))
   {
      /* For each atom in the region                                    */
      for(p=r->start; p!=NULL && p!=r->end; p=FindNextResidue(p))
      {
         /* If it's already in the residue list, just increment the 
            count
         */
         if((l=InResList(PDBResSpec(p), gResList))!=NULL)
         {
            (l->count)++;
         }
         else
         {
            /* It's not already in the residue list, so allocate more
               space and insert into the list
            */
            if(gResList == NULL)
            {
               INIT(gResList, RESLIST);
               sl = gResList;
            }
            else
            {
               ALLOCNEXT(sl, RESLIST);
            }
            if(sl==NULL)
            {
               return(FALSE);
            }
            
            strcpy(sl->res, PDBResSpec(p));
            sl->count = 1;
         }
      }
   }

   return(TRUE);
}


/************************************************************************/
/*>BOOL DoContactAnalysis(FILE *out)
   ---------------------------------
   Does the actual contact analysis. Runs through the loop zones and the
   atoms within these zones. If any pair of atoms is in contact, the
   residue pair is added to the contact list. Finally the contacts
   seen in the current PDB file are displayed and the records of total
   contact data are updated.

   25.09.96 Original   By: ACRM
*/
BOOL DoContactAnalysis(FILE *out)
{
   REGION *r1, *r2;
   PDB    *p,  *q;


   ClearContacts();
   
   /* For each region                                                   */
   for(r1=gLoops; r1!=NULL; NEXT(r1))
   {
      /* For each other region                                          */
      for(r2=r1->next; r2!=NULL; NEXT(r2))
      {
         for(p=r1->start; p!=r1->end; NEXT(p))
         {
            for(q=r2->start; q!=r2->end; NEXT(q))
            {
               if(DISTSQ(p,q) <= gDistSqMax)
               {
                  /* Store the contact (i.e. the pair of residues       */
                  if(!StoreContact(p, q))
                  {
                     ClearContacts();
                     fprintf(stderr,"No memory to store contact data\n");
                     return(FALSE);
                  }
                  /* Store the contacting residues separately           */
                  if(!StoreCRes(p))
                  {
                     ClearContacts();
                     fprintf(stderr,"No memory to store contact data\n");
                     return(FALSE);
                  }
                  if(!StoreCRes(q))
                  {
                     ClearContacts();
                     fprintf(stderr,"No memory to store contact data\n");
                     return(FALSE);
                  }
                  
                  /* One could try to do something clever here to step p 
                     and q on to the end of the residue to save time, but
                     for now we won't bother!
                  */
               }
            }
         }
      }
   }

   DisplayContacts(out, gContacts, (REAL)(-1.0));
   if(!UpdateTotals())
   {
      fprintf(stderr,"No memory to update total contact data\n");
      return(FALSE);
   }
   
   return(TRUE);
}


/************************************************************************/
/*>void ClearContacts(void)
   ------------------------
   Frees up the contact list which is used to store contacts for a single
   protein.

   25.09.96 Original   By: ACRM
   01.10.96 Also frees the list of residues which make contact
*/
void ClearContacts(void)
{
   FREELIST(gContacts, CONTACT);
   gContacts = NULL;

   FREELIST(gCResList, RESLIST);
   gCResList = NULL;
}


/************************************************************************/
/*>BOOL StoreContact(PDB *p, PDB *q)
   ---------------------------------
   Store a contact between two residues into a linked list of contacts.
   The contact is stored in both directions.

   We set the counts to zero here since we are just storing a list of
   what makes contact; we're not using any counts yet

   25.09.96 Original   By: ACRM
*/
BOOL StoreContact(PDB *p, PDB *q)
{
   static CONTACT *c;
   char   resspec1[16],
          resspec2[16];


   /* Create residue specs from the PDB pointers                        */
   strcpy(resspec1, PDBResSpec(p));
   strcpy(resspec2, PDBResSpec(q));

   if(gContacts == NULL)
   {
      INIT(gContacts, CONTACT);
      c = gContacts;
   }
   else
   {
      if(GotContact(resspec1, resspec2, gContacts)!=NULL)
         return(TRUE);
      
      ALLOCNEXT(c, CONTACT);
   }
   if(c==NULL)
   {
      FREELIST(gContacts, CONTACT);
      return(FALSE);
   }

   /* Store the 1->2 contact                                            */
   strcpy(c->res1, resspec1);
   strcpy(c->res2, resspec2);
   c->count = 0;
   
   /* Allocate space for the 2->1 contact                               */
   ALLOCNEXT(c, CONTACT);
   if(c==NULL)
   {
      FREELIST(gContacts, CONTACT);
      return(FALSE);
   }

   /* Store the 1->2 contact                                            */
   strcpy(c->res1, resspec2);
   strcpy(c->res2, resspec1);
   c->count = 0;

   return(TRUE);
}


/************************************************************************/
/*>void DisplayContacts(FILE *out, CONTACT *clist, REAL cutoff)
   ------------------------------------------------------------
   Displays the contact list such that all contacts made from a given
   residue are grouped. The contacts are only displayed if the fraction
   of contacts made is greater than the specified cutoff

   25.09.96 Original   By: ACRM
   01.10.96 Doesn't print percentages if they are 0.00
*/
void DisplayContacts(FILE *out, CONTACT *clist, REAL cutoff)
{
   CONTACT *c, *d;
   BOOL    Finished;
   int     count, count1, count2;
   REAL    perc;

   
   for(c=clist; c!=NULL; NEXT(c))
   {
      c->done = FALSE;
   }

   Finished = FALSE;
   while(!Finished)
   {
      Finished = TRUE;
      
      for(c=clist; c!=NULL; NEXT(c))
      {
         if(!c->done)
         {
            Finished = FALSE;

            for(d=c; d!=NULL; NEXT(d))
            {
               if(!strcmp(c->res1, d->res1))
               {
                  d->done = TRUE;
                  if((count1 = CountRes(d->res1))==0 || 
                     (count2 = CountRes(d->res2))==0)
                  {
                     count = 1;
                  }
                  else
                  {
                     count = MIN(count1, count2);
                  }
                  perc = (REAL)d->count/(REAL)count;

                  if(perc > cutoff)
                  {
                     if(perc != (REAL)0.0)
                     {
                        fprintf(out,"%-5s makes contact with %-5s %.2f\n",
                                d->res1, d->res2, perc);
                     }
                     else
                     {
                        fprintf(out,"%-5s makes contact with %-5s\n",
                                d->res1, d->res2);
                     }
                  }
               }
            }
         }
      }
   }
}


/************************************************************************/
/*>BOOL UpdateTotals(void)
   -----------------------
   Updates the total contacts list with the contacts seen in this
   individual protein. Also updates the total contacting residue lists.

   25.09.96 Original   By: ACRM
   01.10.96 Added contacting residues list
*/
BOOL UpdateTotals(void)
{
   static CONTACT *t;
   CONTACT        *c, *t2;
   static RESLIST *rt;
   RESLIST        *r, *rt2;
   

   for(c=gContacts; c!=NULL; NEXT(c))
   {
      if(gTotal == NULL)
      {
         INIT(gTotal, CONTACT);
         t = gTotal;

         if(t==NULL)
            return(FALSE);
            
         strcpy(t->res1, c->res1);
         strcpy(t->res2, c->res2);
         t->count = 1;
      }
      else
      {
         /* If we don't have it in the totals, create it                */
         if((t2=GotContact(c->res1, c->res2, gTotal))==NULL)
         {
            ALLOCNEXT(t, CONTACT);
            if(t==NULL)
               return(FALSE);
            
            strcpy(t->res1, c->res1);
            strcpy(t->res2, c->res2);
            t->count = 1;
         }
         else
         {
            strcpy(t2->res1, c->res1);
            strcpy(t2->res2, c->res2);
            (t2->count)++;
         }
      }
   }

   for(r=gCResList; r!=NULL; NEXT(r))
   {
      if(gCRTotal == NULL)
      {
         INIT(gCRTotal, RESLIST);
         rt = gCRTotal;

         if(rt==NULL)
            return(FALSE);
            
         strcpy(rt->res, r->res);
         rt->count = 1;
      }
      else
      {
         /* If we don't have it in the totals, create it                */
         if((rt2=InResList(r->res, gCRTotal))==NULL)
         {
            ALLOCNEXT(rt, RESLIST);
            if(rt==NULL)
               return(FALSE);
            
            strcpy(rt->res, r->res);
            rt->count = 1;
         }
         else
         {
            strcpy(rt2->res, r->res);
            (rt2->count)++;
         }
      }
   }

   return(TRUE);
}


/************************************************************************/
/*>CONTACT *GotContact(char *resspec1, char *resspec2, CONTACT *clist)
   -------------------------------------------------------------------
   Looks in a linked list of contacts to see if a particular contact
   is already found there. If found, the pointer into the linked list
   is returned; if not, returns NULL.

   25.09.96 Original   By: ACRM
*/
CONTACT *GotContact(char *resspec1, char *resspec2, CONTACT *clist)
{
   CONTACT *c;

   
   for(c=clist; c!=NULL; NEXT(c))
   {
      if(!strcmp(c->res1, resspec1) &&
         !strcmp(c->res2, resspec2))
      {
         return(c);
      }
   }
   
   return(NULL);
}


/************************************************************************/
/*>RESLIST *InResList(char *rspec, RESLIST *rlist)
   -----------------------------------------------
   Looks to see if a particular residue specification is seen in the list
   of counts of each residue number. If found, returns a pointer to that
   item in the list; returns NULL if not found.

   30.09.96 Original   By: ACRM
   01.10.96 Added rlist as a parameter
            Changed other input to rspec rather than a PDB pointer
*/
RESLIST *InResList(char *rspec, RESLIST *rlist)
{
   RESLIST *l;
   
   for(l=rlist; l!=NULL; NEXT(l))
   {
      if(!strcmp(rspec, l->res))
         return(l);
   }
   return(NULL);
}


/************************************************************************/
/*>int CountRes(char *rspec)
   -------------------------
   Returns the residue count from the linked list of counts for each
   residue specification for a given residue.

   30.09.96 Original   By: ACRM
*/
int CountRes(char *rspec)
{
   RESLIST *l;

   
   for(l=gResList; l!=NULL; NEXT(l))
   {
      if(!strcmp(rspec, l->res))
         return(l->count);
   }
   
   return(0);
}


/************************************************************************/
/*>char *PDBResSpec(PDB *p)
   ------------------------
   Builds a residue spec of the form [c]nnn[i] from a PDB pointer

   30.09.96 Original   By: ACRM
*/
char *PDBResSpec(PDB *p)
{
   static char rspec[16];
   char        chain[16],
               insert[16];


   strcpy(chain,p->chain);
   if(chain[0] == ' ')
      chain[0] = '\0';
   else
      chain[1] = '\0';

   strcpy(insert,p->insert);
   if(insert[0] == ' ')
      insert[0] = '\0';
   else
      insert[1] = '\0';
   
   sprintf(rspec,"%s%d%s",chain,p->resnum,insert);

   return(rspec);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a usage message.

   30.09.96 Original   By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\nLoopContacts V1.0 (c) 1996, Dr. Andrew C.R. Martin, \
UCL\n");
   fprintf(stderr,"Usage: loopcontacts [controlfile [outputfile]]\n");

   fprintf(stderr,"\nLoopContacts examines contacts between residues in \
a set of loops in\n");
   fprintf(stderr,"one or more proteins. If more than one protein is \
analysed, summary\n");
   fprintf(stderr,"data showing the fraction of times a particular \
residue contact pair\n");
   fprintf(stderr,"is found are generated.\n");

   fprintf(stderr,"\nNote that displayed fractions are relative to the \
number of occurrences\n");
   fprintf(stderr,"of a particular residue number, not just the number \
of proteins.\n");
   fprintf(stderr,"Consequently if residue L30C only occurs in half of \
the structures,\n");
   fprintf(stderr,"but always makes a contact when it does occur, the \
fraction would\n");
   fprintf(stderr,"be displayed as 1.0\n");

   fprintf(stderr,"\nThe program is run from a control file which has \
the following commands:\n");
   fprintf(stderr,"LOOP res1 res2        Specifies the range of a loop. \
Repeated for\n");
   fprintf(stderr,"                      each loop to be examined.\n");
   fprintf(stderr,"PDB filename          Specifies a PDB file to be \
analysed.\n");
   fprintf(stderr,"CONTACTS cutoff       Displays contacts which occur \
in at least cutoff\n");
   fprintf(stderr,"                      fraction of the structures.\n");
   fprintf(stderr,"RESIDUES cutoff       Displays residues which make a \
contact in at least\n");
   fprintf(stderr,"                      cutoff fraction of the \
structures.\n");
   fprintf(stderr,"ECHO text...          Displays arbitrary text in the \
output file.\n");
   fprintf(stderr,"EXIT                  Displays summary contact data \
for all the structures\n");
   fprintf(stderr,"                      analysed (as DISPLAY with a 0.0 \
cutoff) and exits\n");
   fprintf(stderr,"                      the program.\n\n");
   fprintf(stderr,"QUIT                  Exits without display\n");
}


/************************************************************************/
/*>BOOL StoreCRes(PDB *p)
   ----------------------
   Store a residues which makes a contact into a linked list of 
   contacting residues.

   We keep a count of how many contacts are made, though this isn't
   used...

   01.10.96 Original   By: ACRM
*/
BOOL StoreCRes(PDB *p)
{
   static RESLIST *c, *r;

   /* See if p is in the contacting residue list                        */
   if((r=InResList(PDBResSpec(p), gCResList))!=NULL)
   {
      /* Yes, increment the count                                       */
      (r->count)++;
   }
   else
   {
      /* No, allocate space and add to list                             */
      if(gCResList == NULL)
      {
         INIT(gCResList, RESLIST);
         c = gCResList;
      }
      else
      {
         ALLOCNEXT(c, RESLIST);
      }
      if(c==NULL)
      {
         FREELIST(gCResList, RESLIST);
         return(FALSE);
      }

      /* Add to list                                                    */
      strcpy(c->res, PDBResSpec(p));
      c->count = 1;
   }
   
   return(TRUE);
}


/************************************************************************/
void DisplayReslist(FILE *out, RESLIST *rlist, REAL cutoff)
{
   RESLIST *r;
   int     count;
   REAL    perc;

   
   for(r=rlist; r!=NULL; NEXT(r))
   {
      if((count = CountRes(r->res))==0)
      {
         count = 1;
      }
      perc = (REAL)r->count/(REAL)count;

      if(perc > cutoff)
      {
         fprintf(out,"%-5s makes a contact %.2f\n",
                 r->res, perc);
      }
   }
}

