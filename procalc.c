/*************************************************************************

   Program:    ProCalc
   File:       procalc.c
   
   Version:    V1.5
   Date:       11.10.93
   Function:   Simple Protein Calculator
   
   Copyright:  (c) SciTech Software 1993
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0372) 275775
   EMail:      UUCP:  cbmehq!cbmuk!scitec!amartin
                      amartin@scitec.adsp.sub.org
               JANET: andrew@uk.ac.ox.biop
               
**************************************************************************

   This program is not in the public domain, but it may be freely copied
   and distributed for no charge providing this header is included.
   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! The code may not be sold commercially without prior permission 
   from the author, although it may be given away free with commercial 
   products, providing it is made clear that this program is free and that 
   the source code is provided with the program.

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Notes:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  11.10.93 Original

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/pdb.h"
#include "bioplib/parse.h"
#include "bioplib/macros.h"
#include "bioplib/angle.h"

/************************************************************************/
/* Defines
*/
#define MAXNUMPARAM   1
#define MAXSTRPARAM   8
#define MAXSTRLEN     80
#define KEY_PDB       0
#define KEY_DISTANCE  1
#define KEY_ANGLE     2
#define KEY_TORSION   3
#define KEY_HELP      4
#define KEY_QUIT      5
#define NCOMM         6

/************************************************************************/
/* Globals
*/
KeyWd gKeyWords[NCOMM];
char  *gStrParam[MAXSTRPARAM];
REAL  gNumParam[MAXNUMPARAM];

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL SetupParser(void);
PDB *CheckCmdLine(int argc, char **argv);
PDB *OpenAndReadPDB(char *file, FILE *Msgfp);
void DoParseLoop(PDB *pdb);
void ShowDistance(PDB *pdb, char *res1, char *at1, char *res2, char *atm2);
void ShowAngle(PDB *pdb, char *res1, char *at1, char *res2, char *atm2,
               char *res3, char *atm3);
void ShowTorsion(PDB *pdb, char *res1, char *at1, char *res2, char *atm2,
                 char *res3, char *atm3, char *res4, char *atm4);
PDB *GetAtom(PDB *pdb, char *res, char *atm);
void ShowHelp(void);
void Usage(void);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program

   10.10.93 Original    By: ACRM
*/
int main(int argc, char **argv)
{
   PDB *pdb = NULL;
   pdb = CheckCmdLine(argc, argv);
   if(SetupParser())
      DoParseLoop(pdb);
   else
      fprintf(stderr,"Unable to initialise command parser\n");
}

/************************************************************************/
/*>BOOL SetupParser(void)
   ----------------------
   Set up the command parser

   10.10.93 Original    By: ACRM
*/
BOOL SetupParser(void)
{
   int i;

   /* Initialise returned string array                                  */
   for(i=0; i<MAXSTRPARAM; i++)
   {
      if((gStrParam[i] = (char *)malloc(MAXSTRLEN * sizeof(char))) == NULL)
      {
         return(FALSE);
      }
   }

   /* Create keywords                                                   */
   MAKEKEY(gKeyWords[KEY_PDB],      "PDB",      STRING, 1);
   MAKEKEY(gKeyWords[KEY_DISTANCE], "DISTANCE", STRING, 4);
   MAKEKEY(gKeyWords[KEY_ANGLE],    "ANGLE",    STRING, 6);
   MAKEKEY(gKeyWords[KEY_TORSION],  "TORSION",  STRING, 8);
   MAKEKEY(gKeyWords[KEY_HELP],     "HELP",     STRING, 0);
   MAKEKEY(gKeyWords[KEY_QUIT],     "QUIT",     STRING, 0);

   /* Check on MAKEKEYs                                                 */
   for(i=0; i<NCOMM; i++)
      if(gKeyWords[i].name == NULL) return(FALSE);

   return(TRUE);
}

/************************************************************************/
/*>PDB *CheckCmdLine(int argc, char **argv)
   ----------------------------------------
   Check the command line. Read a PDB file if specified or give usage
   message.

   10.10.93 Original    By: ACRM
*/
PDB *CheckCmdLine(int argc, char **argv)
{
   PDB *pdb = NULL;

   if(argc==2)         /* Help (-h) or PDB specified                    */
   {
      if(!strncmp(argv[1],"-h",2))
      {
         Usage();
         exit(0);
      }
      else
      {
         pdb = OpenAndReadPDB(argv[1],stderr);
      }
   }
   else if(argc > 2)
   {
      Usage();
      exit(0);
   }
   return(pdb);
}
/************************************************************************/
/*>PDB *OpenAndReadPDB(char *file, FILE *Msgfp)
   --------------------------------------------
   Opens a PDB file specified by name and reads it with ReadPDB. Error
   messages are sent to the specified file. If NULL, no messages will
   be issued. Returns a pointer to the PDB linked list or NULL if failed.

   10.10.93 Original    By: ACRM
   11.10.93 Added legal NULL Msgfp
*/
PDB *OpenAndReadPDB(char *file, FILE *Msgfp)
{
   FILE *fp = NULL;
   PDB  *pdb = NULL;
   int  natoms;

   if((fp=fopen(file,"r"))==NULL)
   {
      if(Msgfp!=NULL) 
         fprintf(Msgfp,"Unable to open file: %s\n",file);
   }
   else
   {
      pdb = ReadPDB(fp, &natoms);
      if(pdb == NULL || natoms == 0)
      {
         if(Msgfp!=NULL) 
            fprintf(Msgfp,"No atoms read from file: %s\n",file);
         if(pdb!=NULL) FREELIST(pdb,PDB);
         pdb = NULL;
      }
      fclose(fp);
   }
   return(pdb);
}
/************************************************************************/
/*>void DoParseLoop(PDB *pdb)
   --------------------------
   The main command parser loop.

   10.10.93 Original    By: ACRM
*/
void DoParseLoop(PDB *pdb)
{
   char buffer[160];
   int key;

   prompt("ProCalc");

   while(fgets(buffer,160,stdin))
   {
      TERMINATE(buffer);
      key = parse(buffer,NCOMM,gKeyWords,gNumParam,gStrParam);
      switch(key)
      {
      case PARSE_ERRC:
         printf("Unknown command: %s\n",buffer);
         break;
      case PARSE_ERRP:
         printf("Error in parameters for command: %s\n",buffer);
         break;
      case KEY_PDB:
         if(pdb != NULL) FREELIST(pdb, PDB);
         pdb = NULL;
         pdb = OpenAndReadPDB(gStrParam[0], stdout);
         break;
      case KEY_DISTANCE:
         ShowDistance(pdb, gStrParam[0], gStrParam[1],
                           gStrParam[2], gStrParam[3]);
         break;
      case KEY_ANGLE:
         ShowAngle(pdb, gStrParam[0], gStrParam[1],
                        gStrParam[2], gStrParam[3],
                        gStrParam[4], gStrParam[5]);
         break;
      case KEY_TORSION:
         ShowTorsion(pdb, gStrParam[0], gStrParam[1],
                        gStrParam[2], gStrParam[3],
                        gStrParam[4], gStrParam[5],
                        gStrParam[6], gStrParam[7]);
         break;
      case KEY_HELP:
         ShowHelp();
         break;
      case KEY_QUIT:
         return;
      }
      prompt("ProCalc");
   }
}
/************************************************************************/
/*>void ShowDistance(PDB *pdb, char *res1, char *atm1, 
                               char *res2, char *atm2)
   ---------------------------------------------------
   Show the distance between 2 atoms

   10.10.93 Original    By: ACRM
*/
void ShowDistance(PDB *pdb, char *res1, char *atm1, 
                            char *res2, char *atm2)
{
   PDB  *at1,
        *at2;
   REAL dist;

   if(pdb==NULL)
   {
      printf("No PDB file loaded\n");
      return;
   }
   at1 = GetAtom(pdb, res1, atm1);
   at2 = GetAtom(pdb, res2, atm2);
   if(at1!=NULL && at2!=NULL)
   {
      dist = DIST(at1, at2);
      WritePDBRecord(stdout,at1);
      WritePDBRecord(stdout,at2);
      printf("Distance = %lf\n",dist);
   }
}
/************************************************************************/
/*>void ShowAngle(PDB *pdb, char *res1, char *atm1, 
                            char *res2, char *atm2,
                            char *res3, char *atm3)
   ------------------------------------------------
   Show the angle between 3 atoms

   10.10.93 Original    By: ACRM
*/
void ShowAngle(PDB *pdb, char *res1, char *atm1, 
                         char *res2, char *atm2,
                         char *res3, char *atm3)
{
   PDB  *at1,
        *at2,
        *at3;
   REAL ang;

   if(pdb==NULL)
   {
      printf("No PDB file loaded\n");
      return;
   }
   at1 = GetAtom(pdb, res1, atm1);
   at2 = GetAtom(pdb, res2, atm2);
   at3 = GetAtom(pdb, res3, atm3);
   if(at1!=NULL && at2!=NULL && at3!=NULL)
   {
      ang = angle(at1->x, at1->y, at1->z,
                  at2->x, at2->y, at2->z,
                  at3->x, at2->y, at3->z);
      WritePDBRecord(stdout,at1);
      WritePDBRecord(stdout,at2);
      WritePDBRecord(stdout,at3);
      printf("Angle = %lf degrees\n",ang*180.0/PI);
   }
}
/************************************************************************/
/*>void ShowTorsion(PDB *pdb, char *res1, char *atm1, 
                              char *res2, char *atm2,
                              char *res3, char *atm3, 
                              char *res4, char *atm4)
   --------------------------------------------------
   Show the torsion angle between 4 atoms

   10.10.93 Original    By: ACRM
*/
void ShowTorsion(PDB *pdb, char *res1, char *atm1, 
                           char *res2, char *atm2,
                           char *res3, char *atm3, 
                           char *res4, char *atm4)
{
   PDB  *at1,
        *at2,
        *at3,
        *at4;
   REAL torsion;

   if(pdb==NULL)
   {
      printf("No PDB file loaded\n");
      return;
   }
   at1 = GetAtom(pdb, res1, atm1);
   at2 = GetAtom(pdb, res2, atm2);
   at3 = GetAtom(pdb, res3, atm3);
   at4 = GetAtom(pdb, res4, atm4);
   if(at1!=NULL && at2!=NULL && at3!=NULL && at4!=NULL)
   {
      torsion = phi(at1->x, at1->y, at1->z,
                    at2->x, at2->y, at2->z,
                    at3->x, at3->y, at3->z,
                    at4->x, at4->y, at4->z);
      WritePDBRecord(stdout,at1);
      WritePDBRecord(stdout,at2);
      WritePDBRecord(stdout,at3);
      WritePDBRecord(stdout,at4);
      printf("Torsion = %lf degrees\n",torsion*180.0/PI);
   }
}
/************************************************************************/
/*>PDB *GetAtom(PDB *pdb, char *res, char *atm)
   --------------------------------------------
   Get a pointer to an atom in a PDB linked list. Takes a residue spec
   of the form [<chain>]resnum[<insert>]. Chain and insert are optional
   and may be in upper or lower case. The atom spec may also be in either
   case.

   10.10.93 Original    By: ACRM
*/
PDB *GetAtom(PDB *pdb, char *res, char *atm)
{
   PDB  *at = NULL;
   int  resnum;
   char chain[8],
        insert[8],
        atspec[8],
        resspec[8];

   strcpy(atspec,atm);
   strcpy(resspec,res);

   UPPER(resspec);
   UPPER(res);
   UPPER(atspec);

   padterm(atspec,4);

   ParseResSpec(resspec, chain, &resnum, insert);

   for(at=pdb; at!=NULL; NEXT(at))
   {
      if(at->resnum    == resnum    &&
         at->chain[0]  == chain[0]  &&
         at->insert[0] == insert[0])
      {
         padterm(at->atnam,4);
         if(!strncmp(at->atnam,atspec,4))
            return(at);
      }
   }
   printf("Residue %s, atom %s not found\n",res,atspec);
   return(NULL);
}
/************************************************************************/
/*>void ShowHelp(void)
   -------------------
   Give brief help on ProCalc

   10.10.93 Original    By: ACRM
*/
void ShowHelp(void)
{
   printf("ProCalc V1.0 11.10.93 - Protein Calculator\n");
   printf("Copyright (c) 1993, Dr. A.C.R. Martin, SciTech Software, \
DKfz\n");
   printf("This program is freely distributable providing no profit is \
made in so doing\n\n");
   printf("PDB <filename>                                              \
Read a PDB file\n");
   printf("DISTANCE <res1> <at1> <res2> <at2>                          \
Calculate distance\n");
   printf("ANGLE <res1> <at1> <res2> <at2> <res3> <at3>                \
Calculate angle\n");
   printf("TORSION <res1> <at1> <res2> <at2> <res3> <at3> <res4> <at4> \
Calculate torsion\n");
   printf("QUIT                                                        \
Exit ProCalc\n");
}
/************************************************************************/
/*>void Usage(void)
   ----------------
   Give command line usage message for ProCalc

   10.10.93 Original    By: ACRM
*/
void Usage(void)
{
   printf("ProCalc V1.0 11.10.93 - Protein Calculator\n");
   printf("Copyright (c) 1993, Dr. A.C.R. Martin, SciTech Software, \
DKfz\n");
   printf("Usage: procalc [<pdbfile>]\n");
   printf("Type 'help' within ProCalc for further information\n\n");
}

