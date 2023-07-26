/*************************************************************************

   Program:    bond
   File:       bond.c
   
   Version:    V1.0
   Date:       19.09.96
   Function:   Interactively calculate new coordinates for an atom given
               the atom to which it bonds, its old position and the new
               bond length.
   
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
#include <math.h>

/************************************************************************/
/* Defines and macros
*/
#define NewPos(scale, a, b) ((a) + (scale)*((b)-(a)))
#define dist(x1, y1, z1, x2, y2, z2)             \
        sqrt((double) (((x1)-(x2))*((x1)-(x2)) + \
                       ((y1)-(y2))*((y1)-(y2)) + \
                       ((z1)-(z2))*((z1)-(z2))))

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
int main(int argc, char **argv)
{
   float ax, ay, az,
         bx, by, bz,
         cx, cy, cz,
         olddist, newdist,
         scale;

   printf("Enter x y z coordinates for first known atom position (A):  ");
   scanf("%f %f %f", &ax, &ay, &az);

   printf("Enter x y z coordinates for second known atom position (B): ");
   scanf("%f %f %f", &bx, &by, &bz);

   olddist = dist(ax, ay, az, bx, by, bz);
   printf("Length of the AB vector is: %.3f\n", olddist); 
   
   printf("Enter required distance along the AB vector: ");
   scanf("%f", &newdist);

   scale = newdist/olddist;

   printf("New x y z coordinates are: %.3f %.3f %.3f\n",
          NewPos(scale, ax, bx),
          NewPos(scale, ay, by),
          NewPos(scale, az, bz));

   return(0);
}


/************************************************************************/
   
