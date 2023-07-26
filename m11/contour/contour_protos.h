/* Prototypes for functions defined in contour.c */
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
