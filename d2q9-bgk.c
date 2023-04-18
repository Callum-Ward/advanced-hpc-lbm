#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "mpi.h"

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

typedef struct 
{
  int start_row;
  int end_row;
} rank_props;

/*
** function prototypes
*/

/* load params, allocate memory, load obstacles & initialise fluid particle densities */
int initialise(const char *restrict paramfile, const char *restrict obstaclefile,
               t_param *restrict params, t_speed **restrict cells_ptr, float **restrict cells_data, t_speed **restrict tmp_cells_ptr, float **restrict tmp_data,
               int **restrict obstacles_ptr, float **restrict av_vels_ptr, rank_props *restrict rank_p, int rank, int size, int *tot_obs);

/*
** The main calculation methods.
** timestep calls, in order, the functions:
** accelerate_flow(), propagate(), rebound() & collision()
*/
void get_rank_sizes(const int rank, const int size,const int rows, rank_props *work);
float timestep(const t_param params, t_speed *restrict cells, t_speed *restrict tmp_cells, int *restrict obstacles, rank_props *rank_p, int rank);
int accelerate_flow(const t_param params, t_speed *restrict cells, int *restrict obstacles, rank_props *rank_p, int rank);
float collision(const t_param params, t_speed *restrict cells, t_speed *restrict tmp_cells, int *restrict obstacles, rank_props *rank_p, int rank);
int write_values(const t_param params,float *restrict cells, int *restrict obstacles, float *restrict av_vels, rank_props *rank_p);

/* finalise, including freeing up allocated memory */
int finalise(const t_param *restrict params, t_speed **restrict cells_ptr, float **restrict cells_data,  t_speed **restrict tmp_cells_ptr, float **restrict tmp_data,
            int **restrict obstacles_ptr, float **restrict av_vels_ptr);

/* Sum all the densities in the grid.
** The total should remain constant from one timestep to the next. */
float total_density(const t_param params, t_speed *restrict cells);

/* compute average velocity */
float av_velocity(const t_param params, float *restrict cells, int *restrict obstacles, rank_props *rank_p);

/* calculate Reynolds number */
float calc_reynolds(const t_param params, float *restrict cells, int *restrict obstacles,rank_props *rank_p);

/* utility functions */
void die(const char *message, const int line, const char *file);
void usage(const char *exe);

/*
** main program:
** initialise, timestep loop, finalise
*/
int main(int argc, char *argv[])
{
  char *paramfile = NULL;    /* name of the input parameter file */
  char *obstaclefile = NULL; /* name of a the input obstacle file */
  t_param params;            /* struct to hold parameter values */
  t_speed *cells = NULL;     /* grid containing fluid densities */
  float *cells_data = NULL;
  t_speed *tmp_cells = NULL; /* scratch space */
  float *tmp_data = NULL;
  int *obstacles = NULL;                                                                  /* grid indicating which cells are blocked */
  float *av_vals = NULL;                                                             /* a record of the av. velocity computed for each timestep */
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

  int tmp_size, tmp_rank, left, right;
  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &tmp_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &tmp_rank);

  const int rank = tmp_rank;
  const int size = tmp_size;
  if (rank==0) printf("Number of ranks %d\n",size);
  int tot_obs;
  left = (rank ==0) ? size-1 : rank-1;
  right = (rank ==size-1) ? 0 : rank+1;
  rank_props rank_p[size];

  initialise(paramfile, obstaclefile, &params, &cells , &cells_data, &tmp_cells, &tmp_data, &obstacles, &av_vals, rank_p, rank, size, &tot_obs);
  
  //printf("my rank: %d start: %d end: %d left: %d right: %d\n", rank, rank_p[rank].start_row, rank_p[rank].end_row, left, right);
  //printf("rank %d init complete\n",rank);
  //printf(" 1 cells->speed0[0]=%f\n",cells->speed0[0]);
  printf("rank %d start %d end %d \n", rank, rank_p[rank].start_row, rank_p[rank].end_row);

  int odd_rank = (rank%2 == 0) ? 0 : 1;
  int even_size = ((size-1)%2 == 0) ? 1 : 0;
  int tag = 77;
  int work_rows = rank_p[rank].end_row - rank_p[rank].start_row;

  int myobs = 0;
  for (size_t i = 0; i < (work_rows+1) * params.nx; i++) {
    if (!obstacles[i]) myobs+=1;
  }

/*   float w0 = params.density * 4.f / 9.f;
  float w1 = params.density / 9.f;
  float w2 = params.density / 36.f;

 for (size_t i = 0; i < (work_rows+3)*params.nx; i++)
  {
    if (cells->speed0[i] != w0 ) printf("speed 0 wrong\n");
    if (cells->speed1[i] != w1 ) printf("speed 1 wrong\n");
    if (cells->speed2[i] != w1 ) printf("speed 2 wrong\n");
    if (cells->speed3[i] != w1 ) printf("speed 3 wrong\n");
    if (cells->speed4[i] != w1 ) printf("speed 4 wrong\n");
    if (cells->speed5[i] != w2 ) printf("speed 5 wrong\n");
    if (cells->speed6[i] != w2 ) printf("speed 6 wrong\n");
    if (cells->speed7[i] != w2 ) printf("speed 7 wrong\n");
    if (cells->speed8[i] != w2 ) printf("speed 8 wrong\n");
  } */

  //printf("rank %d my obs %d\n", rank, myobs);
  //if (rank==0) printf("obs %d\n",tot_obs);

  /* Init time stops here, compute time starts*/
  gettimeofday(&timstr, NULL);
  init_toc = timstr.tv_sec + (timstr.tv_usec / 1000000.0);
  comp_tic = init_toc;

  if (rank == 0 ) printf("total obs = %d\n",tot_obs);

  t_speed *old = NULL;
  for (int tt = 0; tt < params.maxIters; tt++)
  {

    av_vals[tt] = timestep(params, cells, tmp_cells, obstacles, rank_p, rank) / tot_obs;    
    old = cells; // keep pointer to avoid leak
    cells = tmp_cells;
    tmp_cells = old;
    if (tt % 1000 ==0) {
      //printf("rank %d iteration %d av_vals = %f\n",rank, tt, av_vals[tt]);
    }

    
     /*   if (tt == 0 && rank==0) {
      printf("rank %d finished iteration %d av_vals = %f\n", rank, tt, av_vals[tt]);
      printf("speed0[0]=%f\n",cells->speed0[0]);
      printf("data[0]]=%f\n",cells_data[0]);
    } */
    
    MPI_Datatype speed_rows;
    MPI_Type_vector(9, params.nx, (work_rows + 3) * params.nx, MPI_FLOAT, &speed_rows);
    MPI_Type_commit(&speed_rows);

    // exchange halos
    // 1st exchange, even ranks sendrec right
    // 2nd exchange, even ranks sendrec left
    if (size > 0)
    {
      //printf("rank %d halo exchange\n", rank);

      if (odd_rank) //odd sends left
      {
        MPI_Sendrecv(&cells->speed0[params.nx], 1, speed_rows, left, tag, &cells->speed0[0], 1, speed_rows, left, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      } else if (!odd_rank && rank != size-1) { //if last rank number is even don't try and receive from the right yet
        MPI_Sendrecv(&cells->speed0[(work_rows + 1) * params.nx], 1, speed_rows, right, tag, &cells->speed0[(work_rows + 2) * params.nx], 1, speed_rows, right, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }

      if (!odd_rank && !(rank==0 && even_size)) {
        MPI_Sendrecv(&cells->speed0[params.nx], 1, speed_rows, left, tag, &cells->speed0[0], 1, speed_rows, left, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      } else if (odd_rank) {  //if last rank number is even don't try and receive from the right yet
        MPI_Sendrecv(&cells->speed0[(work_rows + 1) * params.nx], 1, speed_rows, right, tag, &cells->speed0[(work_rows + 2) * params.nx], 1, speed_rows, right, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    
      if (even_size) {
        if (rank==0) {
          MPI_Sendrecv(&cells->speed0[params.nx], 1, speed_rows, left, tag, &cells->speed0[0], 1, speed_rows, left, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else if (rank == size - 1)
        {
          MPI_Sendrecv(&cells->speed0[(work_rows + 1) * params.nx], 1, speed_rows, right, tag, &cells->speed0[(work_rows + 2) * params.nx], 1, speed_rows, right, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
      }
    }
/*      if (rank ==0) {
      printf("top halo at tt %d  = %f\n",tt,tmp_cells->speed0[0]);
    } else if (rank == 3) {
      printf("bottom halo at tt %d = %f\n",tt, tmp_cells->speed0[((work_rows+3) * params.nx) -1]);
    } */



#ifdef DEBUG
    printf("==timestep: %d==\n", tt);
    printf("av velocity: %.12E\n", av_vels[tt]);
    printf("tot density: %.12E\n", total_density(params, cells));
#endif
  }
  if(rank==0) printf("make it past compute\n");

  /* Compute time stops here, collate time starts*/
  gettimeofday(&timstr, NULL);
  comp_toc = timstr.tv_sec + (timstr.tv_usec / 1000000.0);
  col_tic = comp_toc;
  
  float *cells_all;
  int *obstacles_all;

  int *rec_cell_count = (int *)malloc(size * sizeof(int));
  int *rec_obs_count = (int *)malloc(size * sizeof(int));
  int *displs = (int *)malloc(size * sizeof(int));
  int *displs_cells = (int *)malloc(size * sizeof(int));
  int obs_sum=0;
  int cells_sum=0;
  for (size_t i = 0; i < size; i++)
  {
    rec_cell_count[i] = ((rank_p[i].end_row - rank_p[i].start_row) + 3) * params.nx * 9;
    rec_obs_count[i] = ((rank_p[i].end_row - rank_p[i].start_row) + 1) * params.nx;
    displs[i] = obs_sum;
    //if (rank == 0) printf("offset obs at %d = %d", i, displs[i]);
    displs_cells[i] = cells_sum;
    //if (rank == 0) printf("offset cells at %d = %d", i, displs_cells[i]);

    obs_sum += rec_obs_count[i];
    cells_sum += rec_cell_count[i];
  }
  //printf("rank %d sending %d obs\n", rank, rec_obs_count[rank]);
  //printf("rank %d sending %d cells\n", rank, rec_cell_count[rank]);

  //if (rank ==0) printf("rank 0 recieving %d obs\n", params.nx * params.ny);
  //if (rank == 0) printf("rank 0 recieving %d cells\n", (params.ny + size * 2) * params.nx * 9);

  if (size > 1)
  {
    if (rank == 0)
    {
      obstacles_all = (int *)_mm_malloc(sizeof(int) * params.ny * params.nx, 64);
      cells_all = (float *)_mm_malloc(sizeof(float) * (params.ny + size * 2) * params.nx * 9, 64);
      MPI_Reduce(MPI_IN_PLACE, av_vals, params.maxIters, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    else
    {
      MPI_Reduce(av_vals, NULL, params.maxIters, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    if (rank == 0) for (size_t i = 0; i < params.maxIters; i += 1000) printf("av_vels[%d]=%f\n", i, av_vals[i]);

    MPI_Gatherv(obstacles, rec_obs_count[rank], MPI_INT, obstacles_all, rec_obs_count, displs, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gatherv(&cells->speed0[0], rec_cell_count[rank], MPI_FLOAT, cells_all, rec_cell_count, displs_cells, MPI_FLOAT, 0, MPI_COMM_WORLD);
  } 
  for (size_t i = 23000; i < 23010; i++)
  {
    //if (rank == 0) printf("cells[%d] = %f\n", i, cells_all[i]);
  } 
  

  /* Total/collate time stops here.*/
  gettimeofday(&timstr, NULL);
  col_toc = timstr.tv_sec + (timstr.tv_usec / 1000000.0);
  tot_toc = col_toc;

  if (rank == 0) {
        /* write final values and free memory */
    printf("==done==\n");
    printf("Reynolds number:\t\t%.12E\n", calc_reynolds(params, cells_all, obstacles_all, rank_p));

    printf("Elapsed Init time:\t\t\t%.6lf (s)\n", init_toc - init_tic);
    printf("Elapsed Compute time:\t\t\t%.6lf (s)\n", comp_toc - comp_tic);
    printf("Elapsed Collate time:\t\t\t%.6lf (s)\n", col_toc - col_tic);
    printf("Elapsed Total time:\t\t\t%.6lf (s)\n", tot_toc - tot_tic);
    write_values(params, cells_all, obstacles_all, av_vals, rank_p);
    finalise(&params,&cells, &tmp_data, &tmp_cells, &tmp_data, &obstacles, &av_vals);


    _mm_free(cells_all);
    cells_all = NULL;
    _mm_free(obstacles_all);
    obstacles_all = NULL;
  }
  MPI_Finalize();
  return EXIT_SUCCESS;
}

void get_rank_sizes(const int rank, const int size, const int tot_rows, rank_props *rank_work) {
  int section_size = tot_rows / size;

  int remainder = tot_rows - (section_size * size);
  //if (rank == 0) printf("rank size remainder %d\n",remainder);
  for (int r = 0; r < size; ++r)
  {
    rank_work[r].start_row = (r * section_size);
    rank_work[r].end_row = ((r + 1) * section_size) - 1;
    //if (rank == 0) printf("my rank: %d start: %d end: %d\n", r, rank_work[r].start_row, rank_work[r].end_row);
  }
  if (remainder)
  {
    int shift = 0;
    for (size_t i = 0; i < size; i++)
    {
      rank_work[i].start_row += shift;
      if (remainder)
      {
        rank_work[i].end_row += (1 + shift);
        shift += 1;
        remainder -= 1;
      }
      else
      {
        rank_work[i].end_row += shift;
      }
      //if (rank == 0) printf(" --- my rank: %d start: %d end: %d\n", i, rank_work[i].start_row, rank_work[i].end_row);
    }
    
  }
   
}

float timestep(const t_param params, t_speed *restrict cells, t_speed *restrict tmp_cells, int *restrict obstacles, rank_props *rank_p, int rank)
{
  if (rank_p[rank].start_row <= (params.ny - 2) && rank_p[rank].end_row >= (params.ny - 2)) {
    //printf("rank %d accelerating\n",rank);
    accelerate_flow(params, cells, obstacles, rank_p, rank);
  }

  return collision(params, cells, tmp_cells, obstacles, rank_p, rank);
}

int accelerate_flow(const t_param params, t_speed *restrict cells, int *restrict obstacles, rank_props *rank_p, int rank)
{
  /* compute weighting factors */
  const float w1 = params.density * params.accel / 9.f;
  const float w2 = params.density * params.accel / 36.f;

  /* modify the 2nd to last row of the grid */
  const int jj = (params.ny - 2) - rank_p[rank].start_row;

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

#pragma vector aligned
#pragma omp simd
  for (int ii = 0; ii < params.nx; ii++)
  {
    /* if the cell is not occupied and
    ** we don't send a negative density */
    if (!obstacles[ii + jj * params.nx] && (cells->speed3[ii + (jj+1) * params.nx] - w1) > 0.f && (cells->speed6[ii + (jj+1) * params.nx] - w2) > 0.f && (cells->speed7[ii + (jj+1) * params.nx] - w2) > 0.f)
    {
      /* increase 'east-side' densities */
      cells->speed1[ii + (jj+1) * params.nx] += w1;
      cells->speed5[ii + (jj+1) * params.nx] += w2;
      cells->speed8[ii + (jj+1) * params.nx] += w2;
      /* decrease 'west-side' densities */
      cells->speed3[ii + (jj+1) * params.nx] -= w1;
      cells->speed6[ii + (jj+1) * params.nx] -= w2;
      cells->speed7[ii + (jj+1) * params.nx] -= w2;
/*       printf("speed0[%d]=%f\n", ii + (jj + 1) * params.nx, cells->speed0[ii + (jj + 1) * params.nx]);
      printf("speed5[%d]=%f\n", ii + (jj + 1) * params.nx, cells->speed5[ii + (jj + 1) * params.nx]);
      printf("speed8[%d]=%f\n", ii + (jj + 1) * params.nx, cells->speed8[ii + (jj + 1) * params.nx]);
      printf("speed3[%d]=%f\n", ii + (jj + 1) * params.nx, cells->speed3[ii + (jj + 1) * params.nx]);
      printf("speed6[%d]=%f\n", ii + (jj + 1) * params.nx, cells->speed6[ii + (jj + 1) * params.nx]);
      printf("speed7[%d]=%f\n", ii + (jj + 1) * params.nx, cells->speed7[ii + (jj + 1) * params.nx]); */
    }
  }

  return EXIT_SUCCESS;
}

const float c_sq = 1.f / 3.f; /* square of speed of sound */
const float w0 = 4.f / 9.f;   /* weighting factor */
const float w1 = 1.f / 9.f;   /* weighting factor */
const float w2 = 1.f / 36.f;  /* weighting factor */

float collision(const t_param params, t_speed *restrict cells, t_speed *restrict tmp_cells, int *restrict obstacles, rank_props *rank_p, int rank)
{

  /* loop over the cells in the grid
  ** NB the collision step is called after
  ** the propagate step and so values of interest
  ** are in the scratch-space grid */

  float tot_u = 0.f;
  // unsigned int tot_cells = 0;

  for (int jj = 1; jj < (rank_p[rank].end_row-rank_p[rank].start_row) +2; jj++)
  {
    const int y_n = jj + 1;
    const int y_s = jj - 1;

    #pragma omp simd
    for (int ii = 0; ii < params.nx; ii++)
    {

      register const int x_e = (ii == params.nx - 1) ? (0) : (ii + 1);
      register const int x_w = (ii == 0) ? (ii + params.nx - 1) : (ii - 1);
      register float speed5 = cells->speed5[x_w + y_s * params.nx]; /* north-east */
      register float speed2 = cells->speed2[ii + y_s * params.nx];  /* north */
      register float speed6 = cells->speed6[x_e + y_s * params.nx]; /* north-west */
      register float speed1 = cells->speed1[x_w + jj * params.nx];  /* east */
      register float speed0 = cells->speed0[ii + jj * params.nx];   /* central cell, no movement */
      register float speed3 = cells->speed3[x_e + jj * params.nx];  /* west */
      register float speed8 = cells->speed8[x_w + y_n * params.nx]; /* south-east */
      register float speed4 = cells->speed4[ii + y_n * params.nx];  /* south */
      register float speed7 = cells->speed7[x_e + y_n * params.nx]; /* south-west */
    /*   if (rank_p[rank].start_row + (jj - 1) == 126)
      {
        printf("speed5=%f\n", speed5);
        printf("speed2=%f\n", speed2);
        printf("speed6=%f\n", speed6);
        printf("speed1=%f\n", speed1);
        printf("speed0=%f\n", speed0);
        printf("speed3=%f\n", speed3);
        printf("speed8=%f\n", speed8);
        printf("speed4=%f\n", speed4);
        printf("speed7=%f\n", speed7);
      }
  */
      /* don't consider occupied cells */
      if (obstacles[ii + (jj-1) * params.nx])
      {
        //if (rank_p[rank].start_row + (jj-1) == 126 ) printf("obstacle at jj:%d ii:%d\n",(jj-1) +rank_p[rank].start_row ,ii);
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
        register const float local_density = speed0 + speed1 + speed2 + speed3 + speed4 + speed5 + speed6 + speed7 + speed8;
        register const float inv_den = 1 / local_density;
        register const float u_x = inv_den * (speed1 + speed5 + speed8 - (speed3 + speed6 + speed7));
        register const float u_y = inv_den * (speed2 + speed5 + speed6 - (speed4 + speed7 + speed8));
        register const float u_sq = u_x * u_x + u_y * u_y;

        //if (rank_p[rank].start_row + (jj - 1) == 126) printf("u_x %f u_y %f at jj:%d ii:%d\n",u_x, u_y, rank_p[rank].start_row + (jj - 1), ii );

        /* equilibrium densities */
        float d_equ;
        register const float sub = u_sq * 1.5f;
        register const float w0d = local_density * w0;
        register const float w1d = w0d * 0.25f;
        register const float w2d = w1d * 0.25f;
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

        tot_u += sqrtf((u_x * u_x) + (u_y * u_y));
        //if (rank_p[rank].start_row + (jj - 1) == 126) printf("tot_u %f at jj:%d ii:%d\n",tot_u, rank_p[rank].start_row + (jj - 1), ii );

      }
    }
    //if (rank_p[rank].start_row + (jj - 1) == 126) printf("tot_t at jj=126 %f\n", tot_u);
  }

  //printf("rank %d collision =  %f\n",rank,tot_u);
  return tot_u;
}

float av_velocity(const t_param params, float *restrict cells, int *restrict obstacles, rank_props *rank_p)
{

  int tot_cells = 0; /* no. of cells used in calculation */
  float tot_u = 0;   /* accumulated magnitudes of velocity for each cell */

/* loop over all non-blocked cells */

  int current_node = 0;
  int work_len = (rank_p[current_node].end_row + 1) * params.nx;
  int work_rows = rank_p[current_node].end_row - rank_p[current_node].start_row;
  int node_offset = 0;
  int jj_offset = 0;
  for (int jj = 0; jj < params.ny; jj++)
  {
#pragma vector aligned
#pragma omp simd
    for (int ii = 0; ii < params.nx; ii++)
    {

      if (jj > rank_p[current_node].end_row)
      {
        jj_offset = jj;
        node_offset += (work_rows + 3) * params.nx * 9;
        current_node += 1;
        work_rows = rank_p[current_node].end_row - rank_p[current_node].start_row;
      }
      int index = ii + (jj - jj_offset) * params.nx;
      /* an occupied cell */
      if (!obstacles[ii + jj * params.nx])
      {
        /* local density total */
        float local_density = 0.f;
        int baseIndex = node_offset + params.nx + index;
        int speedOffset = (work_rows + 3) *params.nx;

        local_density += cells[baseIndex];
        local_density += cells[baseIndex + speedOffset];
        local_density += cells[baseIndex + speedOffset*2];
        local_density += cells[baseIndex + speedOffset*3];
        local_density += cells[baseIndex + speedOffset*4];
        local_density += cells[baseIndex + speedOffset*5];
        local_density += cells[baseIndex + speedOffset*6];
        local_density += cells[baseIndex + speedOffset*7];
        local_density += cells[baseIndex + speedOffset*8];
        float inv = 1 / local_density;

        /* x-component of velocity */
        float u_x = (cells[baseIndex + speedOffset] + cells[baseIndex + speedOffset * 5] + cells[baseIndex + speedOffset * 8] - (cells[baseIndex + speedOffset * 3] + cells[baseIndex + speedOffset * 6] + cells[baseIndex + speedOffset*7])) * inv;
        /* compute y velocity component */
        float u_y = (cells[baseIndex + speedOffset * 2] + cells[baseIndex + speedOffset * 5] + cells[baseIndex + speedOffset * 6] - (cells[baseIndex + speedOffset * 4] + cells[baseIndex + speedOffset * 7] + cells[baseIndex + speedOffset*8])) * inv;
        /* accumulate the norm of x- and y- velocity components */
        tot_u += sqrtf((u_x * u_x) + (u_y * u_y));
        /* increase counter of inspected cells */
        ++tot_cells;
      }
    }
  } 

    return tot_u / (float)tot_cells;
  }

int initialise(const char *restrict paramfile, const char *restrict obstaclefile,
               t_param *params, t_speed **restrict cells_ptr, float **restrict cells_data, t_speed **restrict tmp_cells_ptr, float **restrict tmp_data,
                int **restrict obstacles_ptr, float **restrict av_vels_ptr, rank_props *restrict rank_p, int rank, int size, int *tot_obs)
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
  get_rank_sizes(rank, size, params->ny, rank_p); //populate rank_p struct to with relevent processing regions per rank
  int work_rows = rank_p[rank].end_row - rank_p[rank].start_row + 3;//allocate data for rank processing region & 2 halos
  for (size_t i = 0; i < size; i++)
  {
    //if (rank ==0) printf("rank %d from %d to %d\n",i,rank_p[i].start_row, rank_p[i].end_row);
  }
  
  //printf("my rank %d work %d\n",rank ,work_rows);

  
  /* main grid */
  *cells_data = (float *)_mm_malloc(sizeof(float) * work_rows * params->nx * 9, 64);
  *cells_ptr = (t_speed *)_mm_malloc(sizeof(t_speed) , 64);
  (*cells_ptr)->speed0 = &(*cells_data[0]);
  (*cells_ptr)->speed1 = &(*cells_data)[(work_rows) * params->nx];
  (*cells_ptr)->speed2 = &(*cells_data)[(work_rows) * params->nx*2];
  (*cells_ptr)->speed3 = &(*cells_data)[(work_rows) * params->nx*3];
  (*cells_ptr)->speed4 = &(*cells_data)[(work_rows) * params->nx*4];
  (*cells_ptr)->speed5 = &(*cells_data)[(work_rows) * params->nx*5];
  (*cells_ptr)->speed6 = &(*cells_data)[(work_rows) * params->nx*6];
  (*cells_ptr)->speed7 = &(*cells_data)[(work_rows) * params->nx*7];
  (*cells_ptr)->speed8 = &(*cells_data)[(work_rows) * params->nx*8];

/*   (*cells_ptr)->speed0[params->nx] = 44.f;
  printf("data[params.nx]=%f\n", (*cells_data)[params->nx]);
  (*cells_ptr)->speed0[(work_rows+3) * params->nx + params->nx] = 54.f;
  printf("data[work_rows+3) * params->nx + params.nx]=%f\n", (*cells_data)[(work_rows+3) * params->nx + params->nx]); */

  if (*cells_ptr == NULL)
    die("cannot allocate memory for cells", __LINE__, __FILE__);

  /* 'helper' grid, used as scratch space */
  *tmp_data = (float *)_mm_malloc(sizeof(float) * work_rows * params->nx * 9, 64);
  *tmp_cells_ptr = (t_speed *)_mm_malloc(sizeof(t_speed) , 64);
  (*tmp_cells_ptr)->speed0 = &(*tmp_data)[0];
  (*tmp_cells_ptr)->speed1 = &(*tmp_data)[(work_rows) * params->nx];
  (*tmp_cells_ptr)->speed2 = &(*tmp_data)[(work_rows) * params->nx * 2];
  (*tmp_cells_ptr)->speed3 = &(*tmp_data)[(work_rows) * params->nx * 3];
  (*tmp_cells_ptr)->speed4 = &(*tmp_data)[(work_rows) * params->nx * 4];
  (*tmp_cells_ptr)->speed5 = &(*tmp_data)[(work_rows) * params->nx * 5];
  (*tmp_cells_ptr)->speed6 = &(*tmp_data)[(work_rows) * params->nx * 6];
  (*tmp_cells_ptr)->speed7 = &(*tmp_data)[(work_rows) * params->nx * 7];
  (*tmp_cells_ptr)->speed8 = &(*tmp_data)[(work_rows) * params->nx * 8];

  if (*tmp_cells_ptr == NULL)
    die("cannot allocate memory for tmp_cells", __LINE__, __FILE__);

  /* the map of obstacles */
  *obstacles_ptr = (int *)_mm_malloc(sizeof(int)*(work_rows-2)*params->nx, 64);
  //printf("my rank %d obs %d\n", rank, work_rows - 2);

  if (*obstacles_ptr == NULL)
    die("cannot allocate column memory for obstacles", __LINE__, __FILE__);

  /* initialise densities */
  float w0 = params->density * 4.f / 9.f;
  float w1 = params->density / 9.f;
  float w2 = params->density / 36.f;

#pragma omp simd
  for (int jj = 0; jj < work_rows; jj++) //write to cells including halo regions
  {
    for (int ii = 0; ii < params->nx; ii++)
    {
      //sprintf("rank %d cells_data[%d]=%f\n", rank, ii + jj * params->nx, cells_data[ii + jj * params->nx]);
      /* centre */
      (*cells_ptr)->speed0[ii + jj * params->nx] = w0;
      /* axis directions */
      (*cells_ptr)->speed1[ii + jj * params->nx] = w1;
      (*cells_ptr)->speed2[ii + jj * params->nx] = w1;
      (*cells_ptr)->speed3[ii + jj * params->nx] = w1;
      (*cells_ptr)->speed4[ii + jj * params->nx] = w1;

      /* diagonals */
      (*cells_ptr)->speed5[ii + jj * params->nx] = w2;
      (*cells_ptr)->speed6[ii + jj * params->nx] = w2;
      (*cells_ptr)->speed7[ii + jj * params->nx] = w2;
      (*cells_ptr)->speed8[ii + jj * params->nx] = w2;

       /* first set all cells in obstacle array to zero */
      if (ii + jj * params->nx < (work_rows -2)*params->nx)
      {
        (*obstacles_ptr)[ii + jj * params->nx] = 0;
      }
     
    }
  }
  //printf("cells_ptr->speed0[0]=%f\n",(*cells_ptr)->speed0[0]);

  /* open the obstacle data file */
  fp = fopen(obstaclefile, "r");

  if (fp == NULL)
  {
    sprintf(message, "could not open input obstacles file: %s", obstaclefile);
    die(message, __LINE__, __FILE__);
  }
  /* read-in the blocked cells list */
  *tot_obs = 0;
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

    if (blocked) *tot_obs +=1;
    if (yy >= rank_p[rank].start_row && yy <= rank_p[rank].end_row && blocked) { //only read obstacles relevant to rank processing region
      (*obstacles_ptr)[xx + (yy-rank_p[rank].start_row) * params->nx] = blocked;
 
    }
  }
  /* and close the file */
  *tot_obs = (params->ny*params->nx) - *tot_obs;

  fclose(fp);

  /*
  ** allocate space to hold a record of the avarage velocities computed
  ** at each timestep
  */
  *av_vels_ptr = (float *)_mm_malloc(sizeof(float) * params->maxIters, 64);

  return EXIT_SUCCESS;
}

int finalise(const t_param *params, t_speed **restrict cells_ptr, float **restrict cells_data, t_speed **restrict tmp_cells_ptr, float **restrict tmp_data,
             int **restrict obstacles_ptr, float **av_vels_ptr)
{
  /*
  ** _mm_free up allocated memory
  */

  _mm_free(*cells_ptr);
  *cells_ptr = NULL;

  _mm_free(*cells_data);
  *cells_data = NULL;

  _mm_free(*tmp_cells_ptr);
  *tmp_cells_ptr = NULL;

  _mm_free(*tmp_data);
  *tmp_data = NULL;

  _mm_free(*obstacles_ptr);
  *obstacles_ptr = NULL;

  _mm_free(*av_vels_ptr);
  *av_vels_ptr = NULL;
  return EXIT_SUCCESS;
}

float calc_reynolds(const t_param params, float *restrict cells, int *restrict obstacles, rank_props *rank_p)
{
  const float viscosity = 1.f / 6.f * (2.f / params.omega - 1.f);

  return av_velocity(params, cells, obstacles, rank_p) * params.reynolds_dim / viscosity;
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

int write_values(const t_param params, float *restrict cells, int *restrict obstacles, float *av_vels, rank_props *rank_p)
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

  int current_node = 0;
  int work_len = (rank_p[current_node].end_row + 1) * params.nx;
  int work_rows = rank_p[current_node].end_row - rank_p[current_node].start_row;
  int node_offset=0;
  int jj_offset = 0;
  for (int jj = 0; jj < params.ny; jj++)
  {
    for (int ii = 0; ii < params.nx; ii++)
    {
      
      if (jj > rank_p[current_node].end_row) {
        jj_offset = jj;
        node_offset += (work_rows + 3) * params.nx * 9;
        current_node += 1;
        work_rows = rank_p[current_node].end_row - rank_p[current_node].start_row;
      }
      int index = ii + (jj-jj_offset) * params.nx;
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
        int baseIndex = node_offset + params.nx + index;
        int speedOffset = (work_rows + 3) *params.nx;

        local_density += cells[baseIndex];
        local_density += cells[baseIndex + speedOffset];
        local_density += cells[baseIndex + speedOffset *2];
        local_density += cells[baseIndex + speedOffset *3];
        local_density += cells[baseIndex + speedOffset *4];
        local_density += cells[baseIndex + speedOffset *5];
        local_density += cells[baseIndex + speedOffset *6];
        local_density += cells[baseIndex + speedOffset *7];
        local_density += cells[baseIndex + speedOffset *8];


        /* compute x velocity component */
        // order 1,5,8,3,6,7
        u_x = (cells[baseIndex + speedOffset] + cells[baseIndex + speedOffset * 5] + cells[baseIndex + speedOffset * 8] - (cells[baseIndex + speedOffset * 3] + cells[baseIndex + speedOffset * 6] + cells[baseIndex + speedOffset*7])) / local_density;
        /* compute y velocity component */
        //order 2,5,6,4,7,8
        u_y = (cells[baseIndex + speedOffset * 2] + cells[baseIndex + speedOffset * 5] + cells[baseIndex + speedOffset * 6] - (cells[baseIndex + speedOffset * 4] + cells[baseIndex + speedOffset * 7] + cells[baseIndex + speedOffset*8])) / local_density;
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
