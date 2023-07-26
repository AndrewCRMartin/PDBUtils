#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/MathType.h"


void CalcCellTrans(VEC3F UnitCell, VEC3F CellAngles, 
                   VEC3F *xtrans, VEC3F *ytrans, VEC3F *ztrans);
PDB *DupePDB(PDB *in);
void SetCoordinates(PDB *out, PDB *in);



int main(int argc, char **argv)
{
   PDB   *pdb,
         *pdb2;
   int   natom;
   VEC3F UnitCell,
         CellAngles,
         xtrans, ytrans, ztrans;
   FILE  *fp;

   
   

   fp = fopen("hy5_4a.pdb","r");
   pdb = ReadPDB(fp,&natom);
   
   UnitCell.x = 54.790;
   UnitCell.y = 74.820;
   UnitCell.z = 78.950;
   CellAngles.x = 90.00 * PI / 180.0;
   CellAngles.y = 101.82 * PI / 180.0;
   CellAngles.z = 90.00 * PI / 180.0;
   
   CalcCellTrans(UnitCell, CellAngles, &xtrans, &ytrans, &ztrans);
   
   WritePDB(stdout,pdb);
   
   pdb2 = DupePDB(pdb);
   TranslatePDB(pdb2, xtrans);
   WritePDB(stdout,pdb2);

   SetCoordinates(pdb2, pdb);
   TranslatePDB(pdb2, ytrans);
   WritePDB(stdout,pdb2);
   
   SetCoordinates(pdb2, pdb);
   TranslatePDB(pdb2, ztrans);
   WritePDB(stdout,pdb2);

   xtrans.x *= (-1);
   xtrans.y *= (-1);
   xtrans.z *= (-1);
   ytrans.x *= (-1);
   ytrans.y *= (-1);
   ytrans.z *= (-1);
   ztrans.x *= (-1);
   ztrans.y *= (-1);
   ztrans.z *= (-1);

   SetCoordinates(pdb2, pdb);
   TranslatePDB(pdb2, xtrans);
   WritePDB(stdout,pdb2);

   SetCoordinates(pdb2, pdb);
   TranslatePDB(pdb2, ytrans);
   WritePDB(stdout,pdb2);
   
   SetCoordinates(pdb2, pdb);
   TranslatePDB(pdb2, ztrans);
   WritePDB(stdout,pdb2);

   return(0);
}



void SetCoordinates(PDB *out, PDB *in)
{
   PDB *p, *q;
   
   for(p=in, q=out; p!=NULL && q!=NULL; NEXT(p), NEXT(q))
   {
      q->x = p->x;
      q->y = p->y;
      q->z = p->z;
   }
}

      
   
   



void CalcCellTrans(VEC3F UnitCell, VEC3F CellAngles, 
                   VEC3F *xtrans, VEC3F *ytrans, VEC3F *ztrans)
{
   REAL lena, lenb, lenc,
        cosa, cosb, cosg, sing,
        tmpx, tmpy, tmpz, temp;
   
#ifdef OLD
   lena =  250.0*UnitCell.x; /* A */
   lenb =  250.0*UnitCell.y; /* B */
   lenc = -250.0*UnitCell.z; /* C */
#else
   lena =  UnitCell.x; /* A */
   lenb =  UnitCell.y; /* B */
   lenc = -1.0*UnitCell.z; /* C */
#endif

   cosa = (REAL)cos((double)(CellAngles.x));     /* Alpha */
   cosb = (REAL)cos((double)(CellAngles.y));  /* Beta */
   cosg = (REAL)cos((double)(CellAngles.z));  /* Gamma */
   sing = (REAL)sin((double)(CellAngles.z)); /* Gamma */

   temp = cosa*cosa + cosb*cosb + cosg*cosg - 2.0*cosa*cosb*cosg;
   tmpx = cosb; 
   tmpy = (cosa - cosb*cosg)/sing;
   tmpz = (REAL)sqrt((double)(1.0-temp))/sing;
   
   xtrans->x = lena;
   xtrans->y = (REAL)0.0;
   xtrans->z = (REAL)0.0;
   
   ytrans->x = lenb*cosg;
   ytrans->y = lenb*sing;
   ytrans->z = (REAL)0.0;
   
   ztrans->x = lenc*tmpx;
   ztrans->y = lenc*tmpy;
   ztrans->z = lenc*tmpz;
}

