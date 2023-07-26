/*
   Data contouring from DDJ June 1992, p.91
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"

#if defined (NEVER)
#include <ieeefp.h>
#else
#define NaN 0xFFFFFFFF
#define isnanf(x) ((x) == NaN)
#endif

#define DEFAULT_LEVELS 16

/* Mnemonics for contour line bearings */
#define EAST   0
#define NORTH  1
#define WEST   2
#define SOUTH  3

/* Mnemonics for relative data point positions */
#define SAME      0
#define NEXT      1
#define OPPOSITE  2
#define ADJACENT  3

/* Bit-mapped information in 'map' field        */
#define EW_MAP 0x01
#define NS_MAP 0x02

typedef struct
{
   float x,
         y;
}  LIST;

typedef struct
{
   float    max_value,
            min_value,
            mean,
            std,
            first_level,
            step;
   float    *data;
   char     *map;
   LIST     *list;
   SHORT    dim_x,
            dim_y;
   SHORT    contour_mode;
   USHORT   count;
   char     format[20];
}  GRID;

typedef struct
{
   SHORT    x,
            y;
   UCHAR    bearing;
}  POINT;

#define MXY_to_L(g,x,y) ((USHORT) (y) * (g)->dim_x + (USHORT) (x) + 1)
#define XY_to_L(g,x,y)  ((USHORT) (y) * (g)->dim_x + (USHORT) (x))

extern void ContourText(char *s, float x, float y);
extern void Polyline(int n, LIST *list);


extern int gColourPlot;


/* Contour generation */
void Contour(float *data,
             int dim_x,
             int dim_y,
             double inc);
int scaleData(GRID *grid,
              double inc);
static void startLine(GRID *grid);
static void startEdge(GRID *grid,
                      float level,
                      UCHAR bearing);
void startInterior(GRID *grid,
                   float level);
void drawLine(GRID *grid,
              POINT *point,
              float level);
static void markInUse(GRID *grid,
                      POINT *point,
                      UCHAR face);
static UCHAR faceInUse(GRID *grid,
                       POINT *point,
                       UCHAR face);
static void initPoint(GRID *grid);
static void lastPoint(GRID *grid);
static UCHAR savePoint(GRID *grid,
                       POINT *point,
                       float level);
static float getDataPoint(GRID *grid,
                          POINT *point,
                          UCHAR corner);
void SetGrey(GRID *grid, float level);


/************************************************************************/
/*>void Contour(float *data, int dim_x, int dim_y, double inc)
   -----------------------------------------------------------
   Contour a 2D matrix (data) of dimensions dim_x by dim_y
   inc   > 0 increment to use between contour levels.
         < 0 number of contour levels to generate [abs(inc)].
         = 0 generate default number of contour levels.
   
   09.07.92 Typed in.
*/
void Contour(float   *data, 
             int     dim_x, 
             int     dim_y, 
             double  inc)
{
   GRID  grid;
   
   grid.data  = data;
   grid.dim_x = dim_x;
   grid.dim_y = dim_y;
   
   /* Allocate buffers used to contain contour information */
   if((grid.map = (char *)malloc((dim_x + 1) * dim_y)) == NULL)
   {
      fprintf(stderr,"Contour(): unable to allocate buffer! (%d bytes)\n",
              (dim_x + 1) * dim_y * sizeof(LIST));
      return;
   }
   
   if((grid.list = (LIST *)malloc(2 * dim_x * dim_y * sizeof(LIST))) == NULL)
   {
      fprintf(stderr,"Contour(): unable to allocate buffer! (%d bytes)\n",
              2 * dim_x * dim_y * sizeof(LIST));
      free((char *)grid.map);
      return;
   }
   
   /* Generate contours, if not a uniform field. */
   if(scaleData(&grid, inc))
      startLine(&grid);
      
   /* Release memory */
   free((char *)grid.map);
   free((char *)grid.list);
   
   return;
}

/************************************************************************/
/*>int scaleData(GRID *grid, double inc)
   -------------------------------------
   Determine necessary statistics for contouring data set: global max and
   min, etc. Then initialise items used elsewhere.
   
   09.07.92 Typed in.
*/
int scaleData(GRID   *grid,
              double inc)
{
   USHORT   i;
   float    step,
            level,
            sum,
            sum2,
            count,
            p,
            *u,
            *v,
            r;
   char     *s;
   SHORT    n1,
            n2;
   int      first,
            n;
   
   sum = sum2 = count = 0.0;
   
   first = 1;
   s = grid->map;
   u = grid->data;
   v = u + grid->dim_x * grid->dim_y;
   
   for(i=0; i<grid->dim_x * grid->dim_y; i++, u++, v++, s++)
   {
      r = *u;
      sum += r;
      sum2 += r * r;
      count += 1.0;
      
      if(first)
      {
         grid->max_value = grid->min_value = r;
         first = 0;
      }
      else if(grid->max_value < r)
      {
         grid->max_value = r;
      }
      else if(grid->min_value > r)
      {
         grid->min_value = r;
      }
   }
   
   grid->mean = sum / count;
   if(grid->min_value == grid->max_value)
      return(0);
   
   grid->std = sqrt((sum2 - sum * sum / count) / (count - 1.0));
   if(inc > 0.0)
   {
      /* Use specified increment */
      step = inc;
      n = (int)(grid->max_value - grid->min_value) / step + 1;
      
      while(n>40)
      {
         step *= 2.0;
         n = (int)(grid->max_value - grid->min_value) / step + 1;
      }
   }
   else
   {
      /* Choose [specified|reasonable] number of levels and normalise
         increment to a reasonable value.
      */
      n = (inc == 0.0) ? DEFAULT_LEVELS : (SHORT)fabs(inc);
      
      step = 4.0 * grid->std / (float)n;
      p    = pow(10.0, floor(log10((double)step)));
      step = p * floor((step + p / 2.0) / p);
   }
   
   n1 = (int) floor(log10(fabs(grid->max_value)));
   n2 = -((int) floor(log10(step)));
   
   if(n2 > 0)
      sprintf(grid->format, "%%%d.%df", n1+n2+2, n2);
   else
      sprintf(grid->format,"%%%d.0f", n1+1);
   
   if(grid->max_value * grid->min_value < 0.0)
      level = step * floor(grid->mean / step);        /* odd */
   else
      level = step * floor(grid->min_value / step);
   level -= step * floor((float)(n-1)/2);
   
   /* Back up to include add'l levels, if nesc. */
   while(level - step > grid->min_value)
      level -= step;
   
   grid->first_level = level;
   grid->step = step;
   
   return(1);
}

/************************************************************************/
static void startLine(GRID *grid)
{
   USHORT   idx,
            edge;
   double   level;
   
   for(idx=0,level=grid->first_level; 
       level<grid->max_value; 
       level+=grid->step, idx++)
   {
      /* Clear flags */
      grid->contour_mode = (level >= grid->mean);
      memset(grid->map, 0, grid->dim_x * grid->dim_y);
      
      /* Check edges */
      for(edge=0; edge < 4; edge++)
         startEdge(grid, (float)level, edge);
      /* Check interior points */
      startInterior(grid, (float)level);
   }
}

/************************************************************************/
static void startEdge(GRID *grid, float level, UCHAR bearing)
{
   POINT point1,
         point2;
   float last,
         next;
   SHORT i,
         ds;
   
   switch(point1.bearing = bearing)
   {
   case EAST:
      point1.x = 0;
      point1.y = 0;
      ds = 1;
      break;
   case NORTH:
      point1.x = 0;
      point1.y = grid->dim_y - 2;
      ds = 1;
      break;
   case WEST:
      point1.x = grid->dim_x - 2;
      point1.y = grid->dim_y - 2;
      ds = -1;
      break;
   case SOUTH:
      point1.x = grid->dim_x - 2;
      point1.y = 0;
      ds = -1;
      break;
   }
   
   switch(point1.bearing)  /* Find first point with valid data */
   {
   case EAST:
   case WEST:
      next = getDataPoint(grid, &point1, SAME);
      memcpy((char *)&point2, (char *)&point1, sizeof(POINT));
      point2.x -= ds;
      
      for(i=1; i<grid->dim_y; i++, point1.y = point2.y += ds)
      {
         last = next;
         next = getDataPoint(grid, &point1, NEXT);
         if(last >= level && level > next)
         {
            drawLine(grid, &point1, level);
            memcpy((char *)&point1, (char *)&point2, sizeof(POINT));
            point1.x = point2.x + ds;
         }
      }
      break;
   case SOUTH:
   case NORTH:
      next = getDataPoint(grid, &point1, SAME);
      memcpy((char *)&point2, (char *)&point1, sizeof(POINT));
      point2.y += ds;
      
      for(i=1; i<grid->dim_x; i++, point1.x = point2.x += ds)
      {
         last = next;
         next = getDataPoint(grid, &point1, NEXT);
         if(last >= level && level > next)
         {
            drawLine(grid, &point1, level);
            memcpy((char *)&point1, (char *)&point2, sizeof(POINT));
            point1.y = point2.y - ds;
         }
      }
      break;
   }
}

/************************************************************************/
/*>void startInterior(GRID *grid, float level)
   -------------------------------------------
   For specified contour level, check for properly directed contour line 
   for all interior data points. Do _not_ follow contour lines detected by 
   the startEdge() routine.

   10.07.92 Typed
*/
void startInterior(GRID *grid, float level)
{
   POINT point;
   USHORT x,y;
   float next, last;
   
   for(x=1; x<grid->dim_x-1; x++)
   {
      point.x = x;
      point.y = 0;
      point.bearing = EAST;
      next = getDataPoint(grid, &point, SAME);
      for(y=point.y; y<grid->dim_y; y++, point.y++)
      {
         last = next;
         next = getDataPoint(grid, &point, NEXT);
         if(last >= level && level > next)
         {
            if(!faceInUse(grid, &point, WEST))
            {
               drawLine(grid, &point, level);
               point.x = x;
               point.y = y;
               point.bearing = EAST;
            }
         }
      }
   }
}

/************************************************************************/
/*>void drawLine(GRID *grid, POINT *point, float level)
   ----------------------------------------------------
   Given an initial contour point by startEdge() or startInterior(), follow
   contour line until it encounters an edge or a previously contoured cell.
   
   10.07.92 Started typing
   13.07.92 Finished
*/
void drawLine(GRID *grid, POINT *point, float level)
{
   UCHAR exit_bearing;
   UCHAR adj, opp;
   float fadj, fopp;
   
   initPoint(grid);
   
   for(;;)
   {
      /* Add current point to vector list. If either of these points is
         missing, return immediately (open contour)
      */
      if(!savePoint(grid, point, level))
      {
         SetGrey(grid,level);
         lastPoint(grid);
         return;
      }
      
      /* Has this face of this cell been marked in use? If so, this is a
         closed contour.
      */
      if(faceInUse(grid, point, WEST))
      {
         SetGrey(grid,level);
         lastPoint(grid);
         return;
      }
      
      /* Examine adj and opp corners of cell; find appropriate action
      */
      markInUse(grid, point, WEST);
      
      fadj = getDataPoint(grid, point, ADJACENT);
      fopp = getDataPoint(grid, point, OPPOSITE);
      
      /* If either point missing, return immediately (open contour).
      */
      if(isnanf(fadj) || isnanf(fopp))
      {
         SetGrey(grid,level);
         lastPoint(grid);
         return;
      }
      
      adj = (fadj <= level) ? 2 : 0;
      opp = (fopp >= level) ? 1 : 0;
      
      switch(adj + opp)
      {
      case 0:  /* Exit EAST face */
         markInUse(grid, point, NORTH);
         markInUse(grid, point, SOUTH);
         exit_bearing = EAST;
         break;
      case 1:  /* Exit SOUTH face */
         markInUse(grid, point, NORTH);
         markInUse(grid, point, EAST);
         exit_bearing = SOUTH;
         break;
      case 2:  /* Exit NORTH face */
         markInUse(grid, point, EAST);
         markInUse(grid, point, SOUTH);
         exit_bearing = NORTH;
         break;
      case 3:  /* Exit NORTH or SOUTH face, depending on contour level */
         exit_bearing = (grid->contour_mode) ? NORTH : SOUTH;
         break;
      }
      
      /* Update face number, coordinate of defining corner */
      point->bearing = (point->bearing + exit_bearing) % 4;
      switch(point->bearing)
      {
         case EAST:  point->x++; break;
         case NORTH: point->y--; break;
         case WEST:  point->x--; break;
         case SOUTH: point->y++; break;
      }
   }
}

/************************************************************************/
/*>markInUse(GRID *grid, POINT *point, UCHAR face)
   -----------------------------------------------
   Mark specified cell face as contoured. Needed to stop infinte processing
   of closed contours.
   13.07.92 Typed
*/         
static void markInUse(GRID *grid, POINT *point, UCHAR face)
{
   face = (point->bearing + face) % 4;
   switch(face)
   {
   case NORTH:
   case SOUTH:
      grid->map[MXY_to_L(grid, point->x, point->y + (face == SOUTH ? 1 : 0))]
         |= NS_MAP;
      break;
   case EAST:
   case WEST:
      grid->map[MXY_to_L(grid, point->x + (face == EAST ? 1 : 0), point->y)]
         |= EW_MAP;
      break;
   }
}

/************************************************************************/
/*>faceInUse(GRID *grid, POINT *point, UCHAR face)
   -----------------------------------------------
   See if specified cell gace has been marked as contoured.
   13.07.92 Typed
*/
static UCHAR faceInUse(GRID *grid, POINT *point, UCHAR face)
{
   UCHAR r;
   face = (point->bearing + face) % 4;
   switch(face)
   {
   case NORTH:
   case SOUTH:
      r = grid->map[MXY_to_L(grid, point->x, point->y + (face == SOUTH ? 1 : 0))]
         & NS_MAP;
      break;
   case EAST:
   case WEST:
      r = grid->map[MXY_to_L(grid, point->x + (face == EAST ? 1 : 0), point->y)]
         & EW_MAP;
      break;
   }
   return(r);
}

/************************************************************************/
/*>initPoint(GRID *grid)
   ---------------------
   Init contour point list
   13.07.92 Typed
*/
static void initPoint(GRID *grid)
{
   grid->count = 0;
}

/************************************************************************/
/*>lastPoint(GRID *grid)
   ---------------------
   Generate actual contour line from the contour point list.
   13.07.92 Typed
*/
static void lastPoint(GRID *grid)
{
   if(grid->count) Polyline(grid->count, grid->list);
}

/************************************************************************/
/*>savePoint(GRID *grid, POINT *point, float level)
   ------------------------------------------------
   Add specified point to the contour point list
   13.07.92 Typed
*/
static UCHAR savePoint(GRID *grid, POINT *point, float level)
{
   float last, next;
   float x,y,ds;
   char s[80];
   
   static int cnt = 0;
   
   last = getDataPoint(grid, point, SAME);
   next = getDataPoint(grid, point, NEXT);
   
   /* Are points the same value ? */
   if(last == next)
   {
      fprintf(stderr, "(%2d, %2d, %d)  ", point->x, point->y, point->bearing);
      fprintf(stderr, "%8g  %8g  ", last, next);
      fprintf(stderr, "potential divide-by-zero!\n");
      return(0);
   }
   
   x = (float)point->x;
   y = (float)point->y;
   
   ds = (float)((last-level)/(last-next));
   
   switch(point->bearing)
   {
   case EAST:                 y += ds;       break;
   case NORTH: x += ds;       y += 1.0;      break;
   case WEST:  x += 1.0;      y += 1.0 - ds; break;
   case SOUTH: x += 1.0 - ds;                break;
   }
   
   /* Update to contour point list */
   grid->list[grid->count].x = x / (float)(grid->dim_x - 1);
   grid->list[grid->count].y = y / (float)(grid->dim_y - 1);
   
   /* Add text label to contour line */
   if(!(cnt++ % 11))
   {
      sprintf(s, grid->format, level);
      ContourText(s, grid->list[grid->count].x, grid->list[grid->count].y);
   }
   
   /* Update counter */
   grid->count++;
   
   return(1);
}

/************************************************************************/
/*>getDataPoint(GRID *grid, POINT *point, UCHAR corner)
   ----------------------------------------------------
   Return value of datapoint in specified corner of specified cell. point
   parameter contains address of top-left corner of this cell
*/
static float getDataPoint(GRID *grid, POINT *point, UCHAR corner)
{
   UCHAR dx, dy;
   USHORT offset;
   
   switch((point->bearing + corner) % 4)
   {
   case SAME:     dx = 0; dy = 0; break;
   case NEXT:     dx = 0; dy = 1; break;
   case OPPOSITE: dx = 1; dy = 1; break;
   case ADJACENT: dx = 1; dy = 0; break;
   }
   
   offset = XY_to_L(grid, point->x + dx, point->y + dy);
   if((SHORT)(point->x + dx) >= grid->dim_x ||
      (SHORT)(point->y + dy) >= grid->dim_y ||
      (SHORT)(point->x + dx) < 0            ||
      (SHORT)(point->y + dy) < 0)
      return NaN;
   else
      return(grid->data[offset]);
}
   
void SetGrey(GRID *grid, float level)
{
   float r = (float)0.0, 
         g = (float)0.0, 
         b = (float)0.0,
         l;

   l = level/(grid->max_value - grid->min_value);

   if(gColourPlot)
   {
      r = l;
      b = (float)1.0 - l;
      printf("%f %f %f setrgbcolor\n",r,g,b);
   }
   else
   {
      printf("%f setgray\n",level/(grid->max_value - grid->min_value));
   }
}

