/************************************************************************/
/*>VEC3F *AddCB(REAL nx,  REAL ny,  REAL nz, 
                REAL cax, REAL cay, REAL caz, 
                REAL cx,  REAL cy,  REAL cz)
   ------------------------------------------
   Input:    REAL   nx      N  x-coordinate
             REAL   ny      N  y-coordinate
             REAL   nz      N  z-coordinate
             REAL   cax     CA x-coordinate
             REAL   cay     CA y-coordinate
             REAL   caz     CA z-coordinate
             REAL   cz      C  x-coordinate
             REAL   cz      C  y-coordinate
             REAL   cz      C  z-coordinate
   Returns:  VEC3F  *       Pointer to vector containing CB coordinates
   Preconditions: All coordinates must be valid
   Postconditions: The values in the returned VEC3F pointer are only
                   valid until the next call of this routine. It is
                   a static structure from the routine and must not
                   be free()'d.

   Calculates the position for a CB atom given coordinates of N,CA,C
   06.10.98 Original   By: ACRM
*/
VEC3F *AddCB(REAL nx,  REAL ny,  REAL nz, 
             REAL cax, REAL cay, REAL caz, 
             REAL cx,  REAL cy,  REAL cz)
{
    VEC3F  nca, cca, x, y;
    REAL   ang, sx, sy;
    static VEC3F cb;

    /* Calculate vector from N->CA                                      */
    nca.x = cax - nx;
    nca.y = cay - ny;
    nca.z = caz - nz;

    /* Calculate vector from C->CA                                      */
    cca.x = cax - cx;
    cca.y = cay - cy;
    cca.z = caz - cz;

    /* Find sum and cross product of these 2 vectors                    */
    VecAdd3(&x, nca, cca);
    CrossProd3(&y, nca, cca);

    /* Calculate the coordinates for the CB                             */
    /* ang = acos((double)(-1.0))/2.0 - 
             asin((double)(1.0/sqrt((double)3.0))) 
    */
    ang = HALFPI - asin((double)(1.0/sqrt((double)3.0)));
    sx = CACBBOND * cos((double)ang) / VecLen3(x);
    sy = CACBBOND * sin((double)ang) / VecLen3(y);

    cb.x = cax + (x.x * sx) + (y.x * sy);
    cb.y = cay + (x.y * sx) + (y.y * sy);
    cb.z = caz + (x.z * sx) + (y.z * sy);

    return(&cb);
}

