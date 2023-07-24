#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <sys/time.h>
#include <sys/resource.h>
#include <numeric>
#include <execution>
#include <algorithm>
#include <ranges>
#include "/usr/include/experimental/mdspan"
#include "./algos_include/cartesian_product.hpp"

#define NSPEEDS 9
#define FINALSTATEFILE "final_state.dat"
#define AVVELSFILE "av_vels.dat"

/* struct to hold the parameter values */
struct t_param
{
  int nx;           /* no. of cells in x-direction */
  int ny;           /* no. of cells in y-direction */
  int maxIters;     /* no. of iterations */
  int reynolds_dim; /* dimension for Reynolds number */
  float density;    /* density per link */
  float accel;      /* density redistribution */
  float omega;      /* relaxation parameter */
};

/* struct to hold the 'speed' values */
struct t_speed
{
  float *__restrict__ speed0;
  float *__restrict__ speed1;
  float *__restrict__ speed2;
  float *__restrict__ speed3;
  float *__restrict__ speed4;
  float *__restrict__ speed5;
  float *__restrict__ speed6;
  float *__restrict__ speed7;
  float *__restrict__ speed8;
};

/*
** function prototypes
*/

/* load params, allocate memory, load obstacles & initialise fluid particle densities */
int initialise(const char *__restrict__ paramfile, const char *__restrict__ obstaclefile,
               t_param *__restrict__ params, t_speed **__restrict__ cells_ptr, t_speed **__restrict__ tmp_cells_ptr,
               int **__restrict__ obstacles_ptr, float **__restrict__ av_vels_ptr);

/*
** The main calculation methods.
** timestep calls, in order, the functions:
** accelerate_flow(), propagate(), rebound() & collision()
*/
float timestep(const t_param params, t_speed *__restrict__ cells, t_speed *__restrict__ tmp_cells, int *__restrict__ obstacles);
int accelerate_flow(const t_param params, t_speed *__restrict__ cells, int *__restrict__ obstacles);
float collision(const t_param params, t_speed *__restrict__ cells, t_speed *__restrict__ tmp_cells, int *__restrict__ obstacles);
int write_values(const t_param params, t_speed *__restrict__ cells, int *__restrict__ obstacles, float *__restrict__ av_vels);

/* finalise, including freeing up allocated memory */
int finalise(const t_param *__restrict__ params, t_speed **__restrict__ cells_ptr, t_speed **__restrict__ tmp_cells_ptr,
             int **__restrict__ obstacles_ptr, float **__restrict__ av_vels_ptr);

/* Sum all the densities in the grid.
** The total should remain constant from one timestep to the next. */
float total_density(const t_param params, t_speed *__restrict__ cells);

/* compute average velocity */
float av_velocity(const t_param params, t_speed *__restrict__ cells, int *__restrict__ obstacles);

/* calculate Reynolds number */
float calc_reynolds(const t_param params, t_speed *__restrict__ cells, int *__restrict__ obstacles, int obs_count);

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
  int *obstacles = NULL;                                                          /* grid indicating which cells are blocked */
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

  const int obs_count = std::transform_reduce(std::execution::par, &obstacles[0], &obstacles[params.nx*params.ny], 0, std::plus<int>(), [=](int x) { return x == 0 ? 1 : 0;});

  /* Init time stops here, compute time starts*/
  gettimeofday(&timstr, NULL);
  init_toc = timstr.tv_sec + (timstr.tv_usec / 1000000.0);
  comp_tic = init_toc;
  for (int tt = 0; tt < params.maxIters; tt++)
  {
    av_vels[tt] = timestep(params, cells, tmp_cells, obstacles) / obs_count;
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
  printf("Reynolds number:\t\t%.12E\n", calc_reynolds(params, cells, obstacles, obs_count));
  printf("Elapsed Init time:\t\t\t%.6lf (s)\n", init_toc - init_tic);
  printf("Elapsed Compute time:\t\t\t%.6lf (s)\n", comp_toc - comp_tic);
  printf("Elapsed Collate time:\t\t\t%.6lf (s)\n", col_toc - col_tic);
  printf("Elapsed Total time:\t\t\t%.6lf (s)\n", tot_toc - tot_tic);
  write_values(params, cells, obstacles, av_vels);
  finalise(&params, &cells, &tmp_cells, &obstacles, &av_vels);

  return EXIT_SUCCESS;
}

float timestep(const t_param params, t_speed *__restrict__ cells, t_speed *__restrict__ tmp_cells, int *__restrict__ obstacles)
{
  accelerate_flow(params, cells, obstacles);
  return collision(params, cells, tmp_cells, obstacles);
}

int accelerate_flow(const t_param params, t_speed *__restrict__ cells, int *__restrict__ obstacles)
{
  /* compute weighting factors */
  const float w1 = params.density * params.accel / 9.f;
  const float w2 = params.density * params.accel / 36.f;

  /* modify the 2nd row of the grid */
  const int jj = params.ny - 2;

#pragma align cells->speed0 64;
#pragma align cells->speed1 64;
#pragma align cells->speed2 64;
#pragma align cells->speed3 64;
#pragma align cells->speed4 64;
#pragma align cells->speed5 64;
#pragma align cells->speed6 64;
#pragma align cells->speed7 64;
#pragma align cells->speed8 64;
#pragma align obstacles 64;

  auto ids = std::views::common(std::views::iota(jj * params.nx, (jj+1) * params.nx));
  std::for_each(std::execution::par, ids.begin(), ids.end(), [=](int idx) {
    /* if the cell is not occupied and
    ** we don't send a negative density */
    if (!obstacles[idx] && (cells->speed3[idx] - w1) > 0.f && (cells->speed6[idx] - w2) > 0.f && (cells->speed7[idx] - w2) > 0.f)
    {
      /* increase 'east-side' densities */
      cells->speed1[idx] += w1;
      cells->speed5[idx] += w2;
      cells->speed8[idx] += w2;
      /* decrease 'west-side' densities */
      cells->speed3[idx] -= w1;
      cells->speed6[idx] -= w2;
      cells->speed7[idx] -= w2;
    }
  });

  

  return EXIT_SUCCESS;
}

const float w0 = 4.f / 9.f;   /* weighting factor */

float collision(const t_param params, t_speed *__restrict__ cells, t_speed *__restrict__ tmp_cells, int *__restrict__ obstacles)
{

  /* loop over the cells in the grid
  ** NB the collision step is called after
  ** the propagate step and so values of interest
  ** are in the scratch-space grid */

  #pragma align cells->speed0 64;
  #pragma align cells->speed1 64;
  #pragma align cells->speed2 64;
  #pragma align cells->speed3 64;
  #pragma align cells->speed4 64;
  #pragma align cells->speed5 64;
  #pragma align cells->speed6 64;
  #pragma align cells->speed7 64;
  #pragma align cells->speed8 64;
  #pragma align obstacles 64;

  auto xs = std::views::common(std::views::iota(0, params.nx));
  auto ys = std::views::common(std::views::iota(0, params.ny));
  // Combine xs and ys to create 2D range
  auto ids = std::views::cartesian_product(xs, ys);

  return std::transform_reduce(std::execution::par, ids.begin(), ids.end(), 0.0, std::plus<float>(), [=](auto idx) {

    auto [ii, jj] = idx;
    //ii and jj produced from cartesian product are reversed compared to existing indexing logic (ii = row, jj = col) 
    std::swap(ii, jj);

    const int y_n = (jj == params.ny - 1) ? (0) : jj + 1;
    const int y_s = (jj == 0) ? (jj + params.ny - 1) : (jj - 1);
    const int x_e = (ii == params.nx - 1) ? (0) : (ii + 1);
    const int x_w = (ii == 0) ? (ii + params.nx - 1) : (ii - 1);

    float speed5 = cells->speed5[x_w + y_s * params.nx]; /* north-east */
    float speed2 = cells->speed2[ii + y_s * params.nx];  /* north */
    float speed6 = cells->speed6[x_e + y_s * params.nx]; /* north-west */
    float speed1 = cells->speed1[x_w + jj * params.nx];  /* east */
    float speed0 = cells->speed0[ii + jj * params.nx];   /* central cell, no movement */
    float speed3 = cells->speed3[x_e + jj * params.nx];  /* west */
    float speed8 = cells->speed8[x_w + y_n * params.nx]; /* south-east */
    float speed4 = cells->speed4[ii + y_n * params.nx];  /* south */
    float speed7 = cells->speed7[x_e + y_n * params.nx]; /* south-west */

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
      return 0.f;
    }
    else
    {
      /* compute local density total */
      const float local_density = speed0 + speed1 + speed2 + speed3 + speed4 + speed5 + speed6 + speed7 + speed8;
      const float inv_den = 1 / local_density;
      float u_x = inv_den * (speed1 + speed5 + speed8 - (speed3 + speed6 + speed7));
      const float u_y = inv_den * (speed2 + speed5 + speed6 - (speed4 + speed7 + speed8));
      const float u_sq = u_x * u_x + u_y * u_y;

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

      return sqrtf((u_x * u_x) + (u_y * u_y));
    }
  });

}

float av_velocity(const t_param params, t_speed *__restrict__ cells, int *__restrict__ obstacles)
{
  auto ids = std::views::common(std::views::iota(0,params.nx * params.ny));

  return std::transform_reduce(std::execution::par, ids.begin(), ids.end(), 0.0, std::plus<float>(), [=] (int idx) {

    //ii and jj produced from cartesian product are reversed compared to existing idxing logic (ii = row, jj = col)
    
    /* ignore occupied cells */
    if (!obstacles[idx])
    {
      /* local density total */
      float local_density = 0.f;

      local_density += cells->speed0[idx];
      local_density += cells->speed1[idx];
      local_density += cells->speed2[idx];
      local_density += cells->speed3[idx];
      local_density += cells->speed4[idx];
      local_density += cells->speed5[idx];
      local_density += cells->speed6[idx];
      local_density += cells->speed7[idx];
      local_density += cells->speed8[idx];
      const float inv = 1 / local_density;

      /* x-component of velocity */
      const float u_x = (cells->speed1[idx] + cells->speed5[idx] + cells->speed8[idx] - (cells->speed3[idx] + cells->speed6[idx] + cells->speed7[idx])) * inv;
      /* compute y velocity component */
      const float u_y = (cells->speed2[idx] + cells->speed5[idx] + cells->speed6[idx] - (cells->speed4[idx] + cells->speed7[idx] + cells->speed8[idx])) * inv;
      /* accumulate the norm of x- and y- velocity components */
      return sqrtf((u_x * u_x) + (u_y * u_y));
      /* increase counter of inspected cells */
    }
    else return 0.f;
  });

}

int initialise(const char *__restrict__ paramfile, const char *__restrict__ obstaclefile,
               t_param *params, t_speed **__restrict__ cells_ptr, t_speed **__restrict__ tmp_cells_ptr,
               int **__restrict__ obstacles_ptr, float **__restrict__ av_vels_ptr)
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
  //*cells_ptr = (t_speed *)_mm_malloc( sizeof(t_speed) * NSPEEDS,64);

  *cells_ptr = reinterpret_cast<t_speed*>(std::aligned_alloc(64, sizeof(t_speed) * NSPEEDS));

  (*cells_ptr)->speed0 = reinterpret_cast<float*>(std::aligned_alloc(64, sizeof(float) * params->ny * params->nx));
  (*cells_ptr)->speed1 = reinterpret_cast<float*>(std::aligned_alloc(64, sizeof(float) * params->ny * params->nx));
  (*cells_ptr)->speed2 = reinterpret_cast<float*>(std::aligned_alloc(64, sizeof(float) * params->ny * params->nx));
  (*cells_ptr)->speed3 = reinterpret_cast<float*>(std::aligned_alloc(64, sizeof(float) * params->ny * params->nx));
  (*cells_ptr)->speed4 = reinterpret_cast<float*>(std::aligned_alloc(64, sizeof(float) * params->ny * params->nx));
  (*cells_ptr)->speed5 = reinterpret_cast<float*>(std::aligned_alloc(64, sizeof(float) * params->ny * params->nx));
  (*cells_ptr)->speed6 = reinterpret_cast<float*>(std::aligned_alloc(64, sizeof(float) * params->ny * params->nx));
  (*cells_ptr)->speed7 = reinterpret_cast<float*>(std::aligned_alloc(64, sizeof(float) * params->ny * params->nx));
  (*cells_ptr)->speed8 = reinterpret_cast<float*>(std::aligned_alloc(64, sizeof(float) * params->ny * params->nx));

  if (*cells_ptr == NULL)
    die("cannot allocate memory for cells", __LINE__, __FILE__);

  /* 'helper' grid, used as scratch space */
  //*tmp_cells_ptr = (t_speed *)_mm_malloc(sizeof(t_speed) * NSPEEDS,64);

  *tmp_cells_ptr =  reinterpret_cast<t_speed*>(std::aligned_alloc(64, sizeof(t_speed) * NSPEEDS));
  (*tmp_cells_ptr)->speed0 = reinterpret_cast<float*>(std::aligned_alloc(64, sizeof(float) * params->ny * params->nx));
  (*tmp_cells_ptr)->speed1 = reinterpret_cast<float*>(std::aligned_alloc(64, sizeof(float) * params->ny * params->nx));
  (*tmp_cells_ptr)->speed2 = reinterpret_cast<float*>(std::aligned_alloc(64, sizeof(float) * params->ny * params->nx));
  (*tmp_cells_ptr)->speed3 = reinterpret_cast<float*>(std::aligned_alloc(64, sizeof(float) * params->ny * params->nx));
  (*tmp_cells_ptr)->speed4 = reinterpret_cast<float*>(std::aligned_alloc(64, sizeof(float) * params->ny * params->nx));
  (*tmp_cells_ptr)->speed5 = reinterpret_cast<float*>(std::aligned_alloc(64, sizeof(float) * params->ny * params->nx));
  (*tmp_cells_ptr)->speed6 = reinterpret_cast<float*>(std::aligned_alloc(64, sizeof(float) * params->ny * params->nx));
  (*tmp_cells_ptr)->speed7 = reinterpret_cast<float*>(std::aligned_alloc(64, sizeof(float) * params->ny * params->nx));
  (*tmp_cells_ptr)->speed8 = reinterpret_cast<float*>(std::aligned_alloc(64, sizeof(float) * params->ny * params->nx));

  if (*tmp_cells_ptr == NULL)
    die("cannot allocate memory for tmp_cells", __LINE__, __FILE__);

  /* the map of obstacles */

  //*obstacles_ptr = _mm_malloc( sizeof(int) * (params->ny * params->nx),64);
  *obstacles_ptr = reinterpret_cast<int *>(std::aligned_alloc(64, sizeof(int) * (params->ny * params->nx)));

  if (*obstacles_ptr == NULL)
    die("cannot allocate column memory for obstacles", __LINE__, __FILE__);

  /* initialise densities */
  float w0 = params->density * 4.f / 9.f;
  float w1 = params->density / 9.f;
  float w2 = params->density / 36.f;

  auto ids = std::views::common(std::views::iota(0, params->ny * params->nx));
  std::for_each(std::execution::par, ids.begin(), ids.end(), [=](int idx)
  {
    /* centre */
    (*cells_ptr)->speed0[idx] = w0;
    /* axis directions */
    (*cells_ptr)->speed1[idx] = w1;
    (*cells_ptr)->speed2[idx] = w1;
    (*cells_ptr)->speed3[idx] = w1;
    (*cells_ptr)->speed4[idx] = w1;
    /* diagonals */
    (*cells_ptr)->speed5[idx] = w2;
    (*cells_ptr)->speed6[idx] = w2;
    (*cells_ptr)->speed7[idx] = w2;
    (*cells_ptr)->speed8[idx] = w2;
    (*obstacles_ptr)[idx] = 0; 
  });

  // #pragma omp parallel for
  //   for (int jj = 0; jj < params->ny; jj++) {
  //     for (int ii = 0; ii < params->nx; ii++) {

  //       /* centre */
  //       (*cells_ptr)->speed0[ii + jj * params->nx] = w0;
  //       /* axis directions */
  //       (*cells_ptr)->speed1[ii + jj * params->nx] = w1;
  //       (*cells_ptr)->speed2[ii + jj * params->nx] = w1;
  //       (*cells_ptr)->speed3[ii + jj * params->nx] = w1;
  //       (*cells_ptr)->speed4[ii + jj * params->nx] = w1;

  //       /* diagonals */
  //       (*cells_ptr)->speed5[ii + jj * params->nx] = w2;
  //       (*cells_ptr)->speed6[ii + jj * params->nx] = w2;
  //       (*cells_ptr)->speed7[ii + jj * params->nx] = w2;
  //       (*cells_ptr)->speed8[ii + jj * params->nx] = w2;

  //       /* first set all cells in obstacle array to zero */
  //     }
  //   }

  printf(" \n");

  /* open the obstacle data file */
  fp = fopen(obstaclefile, "r");

  if (fp == NULL){
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
  // *av_vels_ptr = (float *)_mm_malloc(sizeof(float) * params->maxIters, 64);
  *av_vels_ptr = reinterpret_cast<float *>(std::aligned_alloc(64, sizeof(float) * params->maxIters));

  return EXIT_SUCCESS;
}

int finalise(const t_param *params, t_speed **__restrict__ cells_ptr, t_speed **__restrict__ tmp_cells_ptr,
             int **__restrict__ obstacles_ptr, float **__restrict__ av_vels_ptr)
{
  /*
  ** _mm_free up allocated memory
  */

  std::free((*cells_ptr)->speed0);
  std::free((*cells_ptr)->speed1);
  std::free((*cells_ptr)->speed2);
  std::free((*cells_ptr)->speed3);
  std::free((*cells_ptr)->speed4);
  std::free((*cells_ptr)->speed5);
  std::free((*cells_ptr)->speed6);
  std::free((*cells_ptr)->speed7);
  std::free((*cells_ptr)->speed8);

  std::free((*tmp_cells_ptr)->speed0);
  std::free((*tmp_cells_ptr)->speed1);
  std::free((*tmp_cells_ptr)->speed2);
  std::free((*tmp_cells_ptr)->speed3);
  std::free((*tmp_cells_ptr)->speed4);
  std::free((*tmp_cells_ptr)->speed5);
  std::free((*tmp_cells_ptr)->speed6);
  std::free((*tmp_cells_ptr)->speed7);
  std::free((*tmp_cells_ptr)->speed8);

  std::free(*cells_ptr);
  *cells_ptr = NULL;

  std::free(*tmp_cells_ptr);
  *tmp_cells_ptr = NULL;

  std::free(*obstacles_ptr);
  *obstacles_ptr = NULL;

  std::free(*av_vels_ptr);
  *av_vels_ptr = NULL;
  return EXIT_SUCCESS;
}

float calc_reynolds(const t_param params, t_speed *__restrict__ cells, int *__restrict__ obstacles, int obs_count)
{
  const float viscosity = 1.f / 6.f * (2.f / params.omega - 1.f);

  return av_velocity(params, cells, obstacles) * params.reynolds_dim / (viscosity * obs_count);
}

float total_density(const t_param params, t_speed *__restrict__ cells)
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

int write_values(const t_param params, t_speed *__restrict__ cells, int *__restrict__ obstacles, float *av_vels)
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
