/*
** Code to implement a d2q9-bgk lattice boltzmann scheme.
** 'd2' inidates a 2-dimensional grid, and
** 'q9' indicates 9 velocities per grid cell.
** 'bgk' refers to the Bhatnagar-Gross-Krook collision step.
**
** The 'speeds' in each cell are numbered as follows:
**
** 6 2 5
**  \|/
** 3-0-1
**  /|\
** 7 4 8
**
** A 2D grid:
**
**           cols
**       --- --- ---
**      | D | E | F |
** rows  --- --- ---
**      | A | B | C |
**       --- --- ---
**
** 'unwrapped' in row major order to give a 1D array:
**
**  --- --- --- --- --- ---
** | A | B | C | D | E | F |
**  --- --- --- --- --- ---
**
** Grid indicies are:
**
**          ny
**          ^       cols(ii)
**          |  ----- ----- -----
**          | | ... | ... | etc |
**          |  ----- ----- -----
** rows(jj) | | 1,0 | 1,1 | 1,2 |
**          |  ----- ----- -----
**          | | 0,0 | 0,1 | 0,2 |
**          |  ----- ----- -----
**          ----------------------> nx
**
** Note the names of the input parameter and obstacle files
** are passed on the command line, e.g.:
**
**   ./d2q9-bgk input.params obstacles.dat
**
** Be sure to adjust the grid dimensions in the parameter file
** if you choose a different obstacle file.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#define NSPEEDS 9
#define FINALSTATEFILE "final_state.dat"
#define AVVELSFILE "av_vels.dat"

/* struct to hold the parameter values */
typedef struct
{
  int nx;           /* no. of cells in x-direction */
  int ny;           /* no. of cells in y-direction */
  int maxIters;     /* no. of iterations */
  int reynolds_dim; /* dimension for Reynolds number */
  float density;    /* density per link */
  float accel;      /* density redistribution */
  float omega;      /* relaxation parameter */
} t_param;

/* struct to hold the 'speed' values */
typedef struct
{
  float *restrict speed0;
  float *restrict speed1;
  float *restrict speed2;
  float *restrict speed3;
  float *restrict speed4;
  float *restrict speed5;
  float *restrict speed6;
  float *restrict speed7;
  float *restrict speed8;
} t_speed;

/*
** function prototypes
*/

/* load params, allocate memory, load obstacles & initialise fluid particle densities */
int initialise(const char *restrict paramfile, const char *restrict obstaclefile,
               t_param *restrict params, t_speed **restrict cells_ptr, t_speed **restrict tmp_cells_ptr,
               int **restrict obstacles_ptr, float **restrict av_vels_ptr);

/*
** The main calculation methods.
** timestep calls, in order, the functions:
** accelerate_flow(), propagate(), rebound() & collision()
*/
float timestep(const t_param params, t_speed *restrict cells, t_speed *restrict tmp_cells, int *restrict obstacles);
int accelerate_flow(const t_param params, t_speed *restrict cells, int *restrict obstacles);
float collision(const t_param params, t_speed *restrict cells, t_speed *restrict tmp_cells, int *restrict obstacles);
int write_values(const t_param params, t_speed *restrict cells, int *restrict obstacles, float *restrict av_vels);

/* finalise, including freeing up allocated memory */
int finalise(const t_param *restrict params, t_speed **restrict cells_ptr, t_speed **restrict tmp_cells_ptr,
             int **restrict obstacles_ptr, float **restrict av_vels_ptr);

/* Sum all the densities in the grid.
** The total should remain constant from one timestep to the next. */
float total_density(const t_param params, t_speed *restrict cells);

/* compute average velocity */
float av_velocity(const t_param params, t_speed *restrict cells, int *restrict obstacles);

/* calculate Reynolds number */
float calc_reynolds(const t_param params, t_speed *restrict cells, int *restrict obstacles);

/* utility functions */
void die(const char *message, const int line, const char *file);
void usage(const char *exe);

/*
** main program:
** initialise, timestep loop, finalise
*/
int main(int argc, char *argv[])
{
  char *paramfile = NULL;                                                            /* name of the input parameter file */
  char *obstaclefile = NULL;                                                         /* name of a the input obstacle file */
  t_param params;                                                                    /* struct to hold parameter values */
  t_speed *cells = NULL;                                                             /* grid containing fluid densities */
  t_speed *tmp_cells = NULL;                                                         /* scratch space */
  int *obstacles = NULL;                                                             /* grid indicating which cells are blocked */
  float *av_vels = NULL;                                                             /* a record of the av. velocity computed for each timestep */
  struct timeval timstr;                                                             /* structure to hold elapsed time */
  double tot_tic, tot_toc, init_tic, init_toc, comp_tic, comp_toc, col_tic, col_toc; /* floating point numbers to calculate elapsed wallclock time */
  /* parse the command line */
  if (argc != 3)
  {
    usage(argv[0]);
  }
  else
  {
    paramfile = argv[1];
    obstaclefile = argv[2];
  }

  /* Total/init time starts here: initialise our data structures and load values from file */
  gettimeofday(&timstr, NULL);
  tot_tic = timstr.tv_sec + (timstr.tv_usec / 1000000.0);
  init_tic = tot_tic;
  initialise(paramfile, obstaclefile, &params, &cells, &tmp_cells, &obstacles, &av_vels);

  /* Init time stops here, compute time starts*/
  gettimeofday(&timstr, NULL);
  init_toc = timstr.tv_sec + (timstr.tv_usec / 1000000.0);
  comp_tic = init_toc;
  for (int tt = 0; tt < params.maxIters; tt++)
  {
    av_vels[tt] = timestep(params, cells, tmp_cells, obstacles);
    t_speed *old = cells; // keep pointer to avoid leak
    cells = tmp_cells;
    tmp_cells = old;
#ifdef DEBUG
    printf("==timestep: %d==\n", tt);
    printf("av velocity: %.12E\n", av_vels[tt]);
    printf("tot density: %.12E\n", total_density(params, cells));
#endif
  }

  /* Compute time stops here, collate time starts*/
  gettimeofday(&timstr, NULL);
  comp_toc = timstr.tv_sec + (timstr.tv_usec / 1000000.0);
  col_tic = comp_toc;

  // Collate data from ranks here

  /* Total/collate time stops here.*/
  gettimeofday(&timstr, NULL);
  col_toc = timstr.tv_sec + (timstr.tv_usec / 1000000.0);
  tot_toc = col_toc;

  /* write final values and free memory */
  printf("==done==\n");
  printf("Reynolds number:\t\t%.12E\n", calc_reynolds(params, cells, obstacles));
  printf("Elapsed Init time:\t\t\t%.6lf (s)\n", init_toc - init_tic);
  printf("Elapsed Compute time:\t\t\t%.6lf (s)\n", comp_toc - comp_tic);
  printf("Elapsed Collate time:\t\t\t%.6lf (s)\n", col_toc - col_tic);
  printf("Elapsed Total time:\t\t\t%.6lf (s)\n", tot_toc - tot_tic);
  write_values(params, cells, obstacles, av_vels);
  finalise(&params, &cells, &tmp_cells, &obstacles, &av_vels);

  return EXIT_SUCCESS;
}

float timestep(const t_param params, t_speed *restrict cells, t_speed *restrict tmp_cells, int *restrict obstacles)
{

  accelerate_flow(params, cells, obstacles);
  // propagate(params, cells, tmp_cells);
  // rebound(params, cells, tmp_cells, obstacles);

  return collision(params, cells, tmp_cells, obstacles);
}

int accelerate_flow(const t_param params, t_speed *restrict cells, int *restrict obstacles)
{
  /* compute weighting factors */
  float w1 = params.density * params.accel / 9.f;
  float w2 = params.density * params.accel / 36.f;

  /* modify the 2nd row of the grid */
  int jj = params.ny - 2;
  int index = jj * params.nx;
#pragma omp simd
  for (int ii = 0; ii < params.nx; ii++)
  {
    /* if the cell is not occupied and
    ** we don't send a negative density */
    if (!obstacles[index] && (cells->speed3[ii + jj * params.nx] - w1) > 0.f && (cells->speed6[ii + jj * params.nx] - w2) > 0.f && (cells->speed7[ii + jj * params.nx] - w2) > 0.f)
    {
      /* increase 'east-side' densities */
      cells->speed1[ii + jj * params.nx] += w1;
      cells->speed5[ii + jj * params.nx] += w2;
      cells->speed8[ii + jj * params.nx] += w2;
      /* decrease 'west-side' densities */
      cells->speed3[ii + jj * params.nx] -= w1;
      cells->speed6[ii + jj * params.nx] -= w2;
      cells->speed7[ii + jj * params.nx] -= w2;
    }
    index++;
  }

  return EXIT_SUCCESS;
}

const float c_sq = 1.f / 3.f; /* square of speed of sound */
const float w0 = 4.f / 9.f;   /* weighting factor */
const float w1 = 1.f / 9.f;   /* weighting factor */
const float w2 = 1.f / 36.f;  /* weighting factor */

float collision(const t_param params, t_speed *restrict cells, t_speed *restrict tmp_cells, int *restrict obstacles)
{

  /* loop over the cells in the grid
  ** NB the collision step is called after
  ** the propagate step and so values of interest
  ** are in the scratch-space grid */
  float local_density;
  float tot_u = 0.f;
  unsigned int tot_cells = 0;

  __assume_aligned(obstacles, 64);
  __assume_aligned(cells->speed0, 64);
  __assume_aligned(cells->speed1, 64);
  __assume_aligned(cells->speed2, 64);
  __assume_aligned(cells->speed3, 64);
  __assume_aligned(cells->speed4, 64);
  __assume_aligned(cells->speed5, 64);
  __assume_aligned(cells->speed6, 64);
  __assume_aligned(cells->speed7, 64);
  __assume_aligned(cells->speed8, 64);

  __assume_aligned(tmp_cells->speed0, 64);
  __assume_aligned(tmp_cells->speed1, 64);
  __assume_aligned(tmp_cells->speed2, 64);
  __assume_aligned(tmp_cells->speed3, 64);
  __assume_aligned(tmp_cells->speed4, 64);
  __assume_aligned(tmp_cells->speed5, 64);
  __assume_aligned(tmp_cells->speed6, 64);
  __assume_aligned(tmp_cells->speed7, 64);
  __assume_aligned(tmp_cells->speed8, 64);
  __assume((params.nx) % 2 == 0);
  __assume((params.ny) % 2 == 0);

  // reduction(+ : tot_cells, tot_u)
  for (int jj = 0; jj < params.ny; jj++)
  {
#pragma omp simd
    for (int ii = 0; ii < params.nx; ii++)
    {

      int y_n = (jj == params.ny - 1) ? (0) : jj + 1;
      int y_s = (jj == 0) ? (jj + params.ny - 1) : (jj - 1);
      int x_e = (ii == params.nx - 1) ? (0) : (ii + 1);
      int x_w = (ii == 0) ? (ii + params.nx - 1) : (ii - 1);
      __assume((obstacles[jj * params.nx + ii]) < 2);
      __assume(y_n < params.ny);
      __assume(x_e < params.nx);
      __assume(y_s < params.ny);
      __assume(x_w < params.nx);
      register float speed5 = cells->speed5[x_w + y_s * params.nx]; /* north-east */
      register float speed2 = cells->speed2[ii + y_s * params.nx];  /* north */
      register float speed6 = cells->speed6[x_e + y_s * params.nx]; /* north-west */
      register float speed1 = cells->speed1[x_w + jj * params.nx];  /* east */
      register float speed0 = cells->speed0[ii + jj * params.nx];   /* central cell, no movement */
      register float speed3 = cells->speed3[x_e + jj * params.nx];  /* west */
      register float speed8 = cells->speed8[x_w + y_n * params.nx]; /* south-east */
      register float speed4 = cells->speed4[ii + y_n * params.nx];  /* south */
      register float speed7 = cells->speed7[x_e + y_n * params.nx]; /* south-west */

      /* don't consider occupied cells */
      if (obstacles[ii + jj * params.nx])
      {
        tmp_cells->speed1[ii + jj * params.nx] = speed3;
        tmp_cells->speed2[ii + jj * params.nx] = speed4;
        tmp_cells->speed3[ii + jj * params.nx] = speed1;
        tmp_cells->speed4[ii + jj * params.nx] = speed2;
        tmp_cells->speed5[ii + jj * params.nx] = speed7;
        tmp_cells->speed6[ii + jj * params.nx] = speed8;
        tmp_cells->speed7[ii + jj * params.nx] = speed5;
        tmp_cells->speed8[ii + jj * params.nx] = speed6;
      }
      else
      {
        /* compute local density total */
        float local_density = speed0 + speed1 + speed2 + speed3 + speed4 + speed5 + speed6 + speed7 + speed8;
        float inv_den = 1 / local_density;
        float u_x = inv_den * (speed1 + speed5 + speed8 - (speed3 + speed6 + speed7));
        float u_y = inv_den * (speed2 + speed5 + speed6 - (speed4 + speed7 + speed8));
        float u_sq = u_x * u_x + u_y * u_y;

        /* equilibrium densities */
        float d_equ;
        const float sub = u_sq * 1.5f;
        const float w0d = local_density * w0;
        const float w1d = w0d * 0.25f;
        const float w2d = w1d * 0.25f;
        /* zero velocity density: weight w0 */
        d_equ = w0d * (1.f - 1.5f * u_sq);
        speed0 = speed0 + params.omega * (d_equ - speed0);
        /* axis speeds: weight w1 */
        d_equ = w1d * (1.f + 3.f * u_x + 4.5f * u_x * u_x - sub);
        speed1 = speed1 + params.omega * (d_equ - speed1);
        d_equ = w1d * (1.f + 3.f * u_y + 4.5f * u_y * u_y - sub);
        speed2 = speed2 + params.omega * (d_equ - speed2);
        d_equ = w1d * (1.f + 3.f * -u_x + 4.5f * -u_x * -u_x - sub);
        speed3 = speed3 + params.omega * (d_equ - speed3);
        d_equ = w1d * (1.f + 3.f * -u_y + 4.5f * -u_y * -u_y - sub);
        speed4 = speed4 + params.omega * (d_equ - speed4);
        /* diagonal speeds: weight w2 */
        d_equ = w2d * (1.f + 3.f * (u_x + u_y) + 4.5f * (u_x + u_y) * (u_x + u_y) - sub);
        speed5 = speed5 + params.omega * (d_equ - speed5);
        d_equ = w2d * (1.f + 3.f * (-u_x + u_y) + 4.5f * (-u_x + u_y) * (-u_x + u_y) - sub);
        speed6 = speed6 + params.omega * (d_equ - speed6);
        d_equ = w2d * (1.f + 3.f * (-u_x - u_y) + 4.5f * (-u_x - u_y) * (-u_x - u_y) - sub);
        speed7 = speed7 + params.omega * (d_equ - speed7);
        d_equ = w2d * (1.f + 3.f * (u_x - u_y) + 4.5f * (u_x - u_y) * (u_x - u_y) - sub);
        speed8 = speed8 + params.omega * (d_equ - speed8);

        tmp_cells->speed0[ii + jj * params.nx] = speed0;
        tmp_cells->speed1[ii + jj * params.nx] = speed1;
        tmp_cells->speed2[ii + jj * params.nx] = speed2;
        tmp_cells->speed3[ii + jj * params.nx] = speed3;
        tmp_cells->speed4[ii + jj * params.nx] = speed4;
        tmp_cells->speed5[ii + jj * params.nx] = speed5;
        tmp_cells->speed6[ii + jj * params.nx] = speed6;
        tmp_cells->speed7[ii + jj * params.nx] = speed7;
        tmp_cells->speed8[ii + jj * params.nx] = speed8;

        // av_velocity step
        // local_density = speed0 + speed1 + speed2 + speed3 + speed4 + speed5 + speed6 + speed7 + speed8;
        // inv_den = 1 / local_density;
        // u_x = inv_den * (speed1 + speed5 + speed8 - speed3 - speed6 - speed7);
        // u_y = inv_den * (speed2 + speed5 + speed6 - speed4 - speed7 - speed8);
        tot_u += sqrtf((u_x * u_x) + (u_y * u_y));
        tot_cells++;
      }
    }
  }
  return tot_u / (float)tot_cells;
}

float av_velocity(const t_param params, t_speed *restrict cells, int *restrict obstacles)
{

  int tot_cells = 0; /* no. of cells used in calculation */
  float tot_u = 0;   /* accumulated magnitudes of velocity for each cell */

/* loop over all non-blocked cells */
#pragma omp simd
  for (int index = 0; index < params.nx * params.ny; index++)
  {

    /* ignore occupied cells */
    if (!obstacles[index])
    {
      /* local density total */
      float local_density = 0.f;

      local_density += cells->speed0[index];
      local_density += cells->speed1[index];
      local_density += cells->speed2[index];
      local_density += cells->speed3[index];
      local_density += cells->speed4[index];
      local_density += cells->speed5[index];
      local_density += cells->speed6[index];
      local_density += cells->speed7[index];
      local_density += cells->speed8[index];
      float inv = 1 / local_density;

      /* x-component of velocity */
      float u_x = (cells->speed1[index] + cells->speed5[index] + cells->speed8[index] - (cells->speed3[index] + cells->speed6[index] + cells->speed7[index])) * inv;
      /* compute y velocity component */
      float u_y = (cells->speed2[index] + cells->speed5[index] + cells->speed6[index] - (cells->speed4[index] + cells->speed7[index] + cells->speed8[index])) * inv;
      /* accumulate the norm of x- and y- velocity components */
      tot_u += sqrtf((u_x * u_x) + (u_y * u_y));
      /* increase counter of inspected cells */
      ++tot_cells;
    }
  }

  return tot_u / (float)tot_cells;
}

int initialise(const char *restrict paramfile, const char *restrict obstaclefile,
               t_param *params, t_speed **restrict cells_ptr, t_speed **restrict tmp_cells_ptr,
               int **restrict obstacles_ptr, float **restrict av_vels_ptr)
{

  char message[1024]; /* message buffer */
  FILE *fp;           /* file pointer */
  int xx, yy;         /* generic array indices */
  int blocked;        /* indicates whether a cell is blocked by an obstacle */
  int retval;         /* to hold return value for checking */

  /* open the parameter file */
  fp = fopen(paramfile, "r");

  if (fp == NULL)
  {
    sprintf(message, "could not open input parameter file: %s", paramfile);
    die(message, __LINE__, __FILE__);
  }

  /* read in the parameter values */
  retval = fscanf(fp, "%d\n", &(params->nx));

  if (retval != 1)
    die("could not read param file: nx", __LINE__, __FILE__);

  retval = fscanf(fp, "%d\n", &(params->ny));

  if (retval != 1)
    die("could not read param file: ny", __LINE__, __FILE__);

  retval = fscanf(fp, "%d\n", &(params->maxIters));

  if (retval != 1)
    die("could not read param file: maxIters", __LINE__, __FILE__);

  retval = fscanf(fp, "%d\n", &(params->reynolds_dim));

  if (retval != 1)
    die("could not read param file: reynolds_dim", __LINE__, __FILE__);

  retval = fscanf(fp, "%f\n", &(params->density));

  if (retval != 1)
    die("could not read param file: density", __LINE__, __FILE__);

  retval = fscanf(fp, "%f\n", &(params->accel));

  if (retval != 1)
    die("could not read param file: accel", __LINE__, __FILE__);

  retval = fscanf(fp, "%f\n", &(params->omega));

  if (retval != 1)
    die("could not read param file: omega", __LINE__, __FILE__);

  /* and close up the file */
  fclose(fp);

  /*
  ** Allocate memory.
  **
  ** Remember C is pass-by-value, so we need to
  ** pass pointers into the initialise function.
  **
  ** NB we are allocating a 1D array, so that the
  ** memory will be contiguous.  We still want to
  ** index this memory as if it were a (row major
  ** ordered) 2D array, however.  We will perform
  ** some arithmetic using the row and column
  ** coordinates, inside the square brackets, when
  ** we want to access elements of this array.
  **
  ** Note also that we are using a structure to
  ** hold an array of 'speeds'.  We will allocate
  ** a 1D array of these structs.
  */

  /* main grid */
  *cells_ptr = (t_speed *)aligned_alloc(64, sizeof(t_speed) * NSPEEDS);

  (*cells_ptr)->speed0 = (float *)aligned_alloc(64, sizeof(float) * params->ny * params->nx);
  (*cells_ptr)->speed1 = (float *)aligned_alloc(64, sizeof(float) * params->ny * params->nx);
  (*cells_ptr)->speed2 = (float *)aligned_alloc(64, sizeof(float) * params->ny * params->nx);
  (*cells_ptr)->speed3 = (float *)aligned_alloc(64, sizeof(float) * params->ny * params->nx);
  (*cells_ptr)->speed4 = (float *)aligned_alloc(64, sizeof(float) * params->ny * params->nx);
  (*cells_ptr)->speed5 = (float *)aligned_alloc(64, sizeof(float) * params->ny * params->nx);
  (*cells_ptr)->speed6 = (float *)aligned_alloc(64, sizeof(float) * params->ny * params->nx);
  (*cells_ptr)->speed7 = (float *)aligned_alloc(64, sizeof(float) * params->ny * params->nx);
  (*cells_ptr)->speed8 = (float *)aligned_alloc(64, sizeof(float) * params->ny * params->nx);

  if (*cells_ptr == NULL)
    die("cannot allocate memory for cells", __LINE__, __FILE__);

  /* 'helper' grid, used as scratch space */
  *tmp_cells_ptr = (t_speed *)aligned_alloc(64, sizeof(t_speed) * NSPEEDS);
  (*tmp_cells_ptr)->speed0 = (float *)aligned_alloc(64, sizeof(float) * params->ny * params->nx);
  (*tmp_cells_ptr)->speed1 = (float *)aligned_alloc(64, sizeof(float) * params->ny * params->nx);
  (*tmp_cells_ptr)->speed2 = (float *)aligned_alloc(64, sizeof(float) * params->ny * params->nx);
  (*tmp_cells_ptr)->speed3 = (float *)aligned_alloc(64, sizeof(float) * params->ny * params->nx);
  (*tmp_cells_ptr)->speed4 = (float *)aligned_alloc(64, sizeof(float) * params->ny * params->nx);
  (*tmp_cells_ptr)->speed5 = (float *)aligned_alloc(64, sizeof(float) * params->ny * params->nx);
  (*tmp_cells_ptr)->speed6 = (float *)aligned_alloc(64, sizeof(float) * params->ny * params->nx);
  (*tmp_cells_ptr)->speed7 = (float *)aligned_alloc(64, sizeof(float) * params->ny * params->nx);
  (*tmp_cells_ptr)->speed8 = (float *)aligned_alloc(64, sizeof(float) * params->ny * params->nx);

  if (*tmp_cells_ptr == NULL)
    die("cannot allocate memory for tmp_cells", __LINE__, __FILE__);

  /* the map of obstacles */
  *obstacles_ptr = aligned_alloc(64, sizeof(int) * (params->ny * params->nx));

  if (*obstacles_ptr == NULL)
    die("cannot allocate column memory for obstacles", __LINE__, __FILE__);

  /* initialise densities */
  float w0 = params->density * 4.f / 9.f;
  float w1 = params->density / 9.f;
  float w2 = params->density / 36.f;
  for (int index = 0; index < params->nx * params->ny; index++)
  {
    /* centre */
    (*cells_ptr)->speed0[index] = w0;
    /* axis directions */
    (*cells_ptr)->speed1[index] = w1;
    (*cells_ptr)->speed2[index] = w1;
    (*cells_ptr)->speed3[index] = w1;
    (*cells_ptr)->speed4[index] = w1;

    /* diagonals */
    (*cells_ptr)->speed5[index] = w2;
    (*cells_ptr)->speed6[index] = w2;
    (*cells_ptr)->speed7[index] = w2;
    (*cells_ptr)->speed8[index] = w2;

    /* first set all cells in obstacle array to zero */
  }
  // #pragma omp simd
  // for (int index = 0; index < params->nx * params->ny; index++) {
  //   (*obstacles_ptr)[index] = 0;
  // }

  /* open the obstacle data file */
  fp = fopen(obstaclefile, "r");

  if (fp == NULL)
  {
    sprintf(message, "could not open input obstacles file: %s", obstaclefile);
    die(message, __LINE__, __FILE__);
  }

  /* read-in the blocked cells list */
  while ((retval = fscanf(fp, "%d %d %d\n", &xx, &yy, &blocked)) != EOF)
  {
    /* some checks */
    if (retval != 3)
      die("expected 3 values per line in obstacle file", __LINE__, __FILE__);

    if (xx < 0 || xx > params->nx - 1)
      die("obstacle x-coord out of range", __LINE__, __FILE__);

    if (yy < 0 || yy > params->ny - 1)
      die("obstacle y-coord out of range", __LINE__, __FILE__);

    if (blocked != 1)
      die("obstacle blocked value should be 1", __LINE__, __FILE__);

    /* assign to array */
    (*obstacles_ptr)[xx + yy * params->nx] = blocked;
  }

  /* and close the file */
  fclose(fp);

  /*
  ** allocate space to hold a record of the avarage velocities computed
  ** at each timestep
  */
  *av_vels_ptr = (float *)malloc(sizeof(float) * params->maxIters);

  return EXIT_SUCCESS;
}

int finalise(const t_param *params, t_speed **restrict cells_ptr, t_speed **restrict tmp_cells_ptr,
             int **restrict obstacles_ptr, float **av_vels_ptr)
{
  /*
  ** free up allocated memory
  */

  free((*cells_ptr)->speed0);
  free((*cells_ptr)->speed1);
  free((*cells_ptr)->speed2);
  free((*cells_ptr)->speed3);
  free((*cells_ptr)->speed4);
  free((*cells_ptr)->speed5);
  free((*cells_ptr)->speed6);
  free((*cells_ptr)->speed7);
  free((*cells_ptr)->speed8);

  free((*tmp_cells_ptr)->speed0);
  free((*tmp_cells_ptr)->speed1);
  free((*tmp_cells_ptr)->speed2);
  free((*tmp_cells_ptr)->speed3);
  free((*tmp_cells_ptr)->speed4);
  free((*tmp_cells_ptr)->speed5);
  free((*tmp_cells_ptr)->speed6);
  free((*tmp_cells_ptr)->speed7);
  free((*tmp_cells_ptr)->speed8);

  free(*cells_ptr);
  *cells_ptr = NULL;

  free(*tmp_cells_ptr);
  *tmp_cells_ptr = NULL;

  free(*obstacles_ptr);
  *obstacles_ptr = NULL;

  free(*av_vels_ptr);
  *av_vels_ptr = NULL;

  return EXIT_SUCCESS;
}

float calc_reynolds(const t_param params, t_speed *restrict cells, int *restrict obstacles)
{
  const float viscosity = 1.f / 6.f * (2.f / params.omega - 1.f);

  return av_velocity(params, cells, obstacles) * params.reynolds_dim / viscosity;
}

float total_density(const t_param params, t_speed *restrict cells)
{
  float total = 0.f; /* accumulator */

  for (int jj = 0; jj < params.ny; jj++)
  {
    for (int ii = 0; ii < params.nx; ii++)
    {
      total += cells->speed0[ii + jj * params.nx];
      total += cells->speed1[ii + jj * params.nx];
      total += cells->speed2[ii + jj * params.nx];
      total += cells->speed3[ii + jj * params.nx];
      total += cells->speed4[ii + jj * params.nx];
      total += cells->speed5[ii + jj * params.nx];
      total += cells->speed6[ii + jj * params.nx];
      total += cells->speed7[ii + jj * params.nx];
      total += cells->speed8[ii + jj * params.nx];
    }
  }

  return total;
}

int write_values(const t_param params, t_speed *restrict cells, int *restrict obstacles, float *av_vels)
{
  FILE *fp;                     /* file pointer */
  const float c_sq = 1.f / 3.f; /* sq. of speed of sound */
  float local_density;          /* per grid cell sum of densities */
  float pressure;               /* fluid pressure in grid cell */
  float u_x;                    /* x-component of velocity in grid cell */
  float u_y;                    /* y-component of velocity in grid cell */
  float u;                      /* norm--root of summed squares--of u_x and u_y */

  fp = fopen(FINALSTATEFILE, "w");

  if (fp == NULL)
  {
    die("could not open file output file", __LINE__, __FILE__);
  }

  for (int jj = 0; jj < params.ny; jj++)
  {
    for (int ii = 0; ii < params.nx; ii++)
    {
      /* an occupied cell */
      if (obstacles[ii + jj * params.nx])
      {
        u_x = u_y = u = 0.f;
        pressure = params.density * c_sq;
      }
      /* no obstacle */
      else
      {
        local_density = 0.f;

        local_density += cells->speed0[ii + jj * params.nx];
        local_density += cells->speed1[ii + jj * params.nx];
        local_density += cells->speed2[ii + jj * params.nx];
        local_density += cells->speed3[ii + jj * params.nx];
        local_density += cells->speed4[ii + jj * params.nx];
        local_density += cells->speed5[ii + jj * params.nx];
        local_density += cells->speed6[ii + jj * params.nx];
        local_density += cells->speed7[ii + jj * params.nx];
        local_density += cells->speed8[ii + jj * params.nx];

        /* compute x velocity component */
        u_x = (cells->speed1[ii + jj * params.nx] + cells->speed5[ii + jj * params.nx] + cells->speed8[ii + jj * params.nx] - (cells->speed3[ii + jj * params.nx] + cells->speed6[ii + jj * params.nx] + cells->speed7[ii + jj * params.nx])) / local_density;
        /* compute y velocity component */
        u_y = (cells->speed2[ii + jj * params.nx] + cells->speed5[ii + jj * params.nx] + cells->speed6[ii + jj * params.nx] - (cells->speed4[ii + jj * params.nx] + cells->speed7[ii + jj * params.nx] + cells->speed8[ii + jj * params.nx])) / local_density;
        /* compute norm of velocity */
        u = sqrtf((u_x * u_x) + (u_y * u_y));
        /* compute pressure */
        pressure = local_density * c_sq;
      }

      /* write to file */
      fprintf(fp, "%d %d %.12E %.12E %.12E %.12E %d\n", ii, jj, u_x, u_y, u, pressure, obstacles[ii + params.nx * jj]);
    }
  }

  fclose(fp);

  fp = fopen(AVVELSFILE, "w");

  if (fp == NULL)
  {
    die("could not open file output file", __LINE__, __FILE__);
  }

  for (int ii = 0; ii < params.maxIters; ii++)
  {
    fprintf(fp, "%d:\t%.12E\n", ii, av_vels[ii]);
  }

  fclose(fp);

  return EXIT_SUCCESS;
}

void die(const char *message, const int line, const char *file)
{
  fprintf(stderr, "Error at line %d of file %s:\n", line, file);
  fprintf(stderr, "%s\n", message);
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void usage(const char *exe)
{
  fprintf(stderr, "Usage: %s <paramfile> <obstaclefile>\n", exe);
  exit(EXIT_FAILURE);
}
