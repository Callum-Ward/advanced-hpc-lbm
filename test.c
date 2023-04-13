#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <mpi.h>

typedef struct
{
  int start_row;
  int end_row;
} work_slice;

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

void get_rank_sizes(const int size, const int tot_rows, work_slice *rank_work)
{
  int section_size = tot_rows / size;
  int remainder = tot_rows - (section_size * size);
  for (int rank = 0; rank < size; ++rank)
  {
    rank_work[rank].start_row = (rank * section_size);
    rank_work[rank].end_row = ((rank + 1) * section_size) - 1;
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
        remainder-=1;
      }
      else
      {
        rank_work[i].end_row += shift;
      }
    }
  }

}

int main(int argc, char *argv[])
{
  int flag, left, right;
  int tmpRank;
  int size;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &tmpRank);

  const int rank = tmpRank;
  int tag = 77;
  left = (rank == 0) ? size - 1 : rank - 1;
  right = (rank == size - 1) ? 0 : rank + 1;
  int odd_rank = (rank % 2 == 0) ? 0 : 1;
  int even_size = ((size - 1) % 2 == 0) ? 1 : 0;
  work_slice work[size];
  int columns = 10;
  int rows = 4;
  get_rank_sizes(size, rows, work);

  int work_rows = work[rank].end_row - work[rank].start_row;


  t_speed *cells = (t_speed *)_mm_malloc(sizeof(t_speed), 64);
  float *data = (float *)_mm_malloc(sizeof(float) * (work_rows+3) * columns * 9, 64);
  // cells->speed0 = (float *)_mm_malloc(sizeof(float *), 64);
  cells->speed0 = &data[0];
  // cells->speed1 = (float *)_mm_malloc(sizeof(float *), 64);
  cells->speed1 = &data[(work_rows+3) * columns];
  // cells->speed2 = (float *)_mm_malloc(sizeof(float *), 64);
  cells->speed2 = &data[(work_rows+3) * columns * 2];
  // cells->speed3 = (float *)_mm_malloc(sizeof(float *), 64);
  cells->speed3 = &data[(work_rows+3) * columns * 3];
  // cells->speed4 = (float *)_mm_malloc(sizeof(float *), 64);
  cells->speed4 = &data[(work_rows+3) * columns * 4];
  // cells->speed5 = (float *)_mm_malloc(sizeof(float *), 64);
  cells->speed5 = &data[(work_rows+3) * columns * 5];
  // cells->speed6 = (float *)_mm_malloc(sizeof(float *), 64);
  cells->speed6 = &data[(work_rows+3) * columns * 6];
  // cells->speed7 = (float *)_mm_malloc(sizeof(float *), 64);
  cells->speed7 = &data[(work_rows+3) * columns * 7];
  // cells->speed8 = (float *)_mm_malloc(sizeof(float *), 64);
  cells->speed8 = &data[(work_rows+3) * columns * 8]; 
  //createCells(cells,cells_data,work_rows+3,10);
  
  float *start = &data[0];
  cells->speed0[0] = 1.f;
  //printf("start[0] = %f data[0] = %f\n",cells->speed0[0],data[0]);

  // for (int i = 0; i < 9* 10 * (work_rows + 3); i++)
  // {
  //     // data[i] = (float)(rank + (i/10*100));
  //     data[i] = (float)rank;
  //     cells->speed0[i+10]
  // }

  for (int i = 0; i < (work_rows+3)  * columns; i++) {
    if (i >= columns && i < (work_rows+2)* columns) {
      cells->speed0[i] = (float)(rank);
      cells->speed1[i] = (float)(rank);
      cells->speed2[i] = (float)(rank);
      cells->speed3[i] = (float)(rank);
      cells->speed4[i] = (float)(rank);
      cells->speed5[i] = (float)(rank);
      cells->speed6[i] = (float)(rank);
      cells->speed7[i] = (float)(rank);
      cells->speed8[i] = (float)(rank);
    } else {
      cells->speed0[i] = (float)0.1*rank ;
      cells->speed1[i] = (float)0.1*rank;
      cells->speed2[i] = (float)0.1*rank;
      cells->speed3[i] = (float)0.1*rank;
      cells->speed4[i] = (float)0.1*rank;
      cells->speed5[i] = (float)0.1*rank;
      cells->speed6[i] = (float)0.1*rank;
      cells->speed7[i] = (float)0.1*rank;
      cells->speed8[i] = (float)0.1*rank;
    }

  


    // cells->speedq[i] = (float)rank+1;
    // cells->speed1[i] = (float)rank+1;
    // cells->speed2[i] = (float)rank+1;
    // cells->speed3[i] = (float)rank+1;
    // cells->speed4[i] = (float)rank+1;
    // cells->speed5[i] = (float)rank+1;
    // cells->speed6[i] = (float)rank+1;
    // cells->speed7[i] = (float)rank+1;
    // cells->speed8[i] = (float)rank+1;
  }

  MPI_Datatype speed_rows;
  MPI_Type_vector(9, columns, (work_rows+3) * columns, MPI_FLOAT, &speed_rows);
  MPI_Type_commit(&speed_rows);
  if (rank==1) {
     for (size_t i = 0; i < 20; i++)
    {
      //printf("cells->speed0[%d]=%f\n",i,cells->speed0[i]);
    }
  }
 
  

  //printf("vector type created\n");
  
  for (int i = 0; i < 1; i++)
  {

    if (odd_rank) // odds send left
    {
      MPI_Sendrecv(&data[columns], 1, speed_rows, left, tag, data, 1, speed_rows, left, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      // printf("rank %d 3\n", rank);
      // printf("round 1 rank= %d send to %d\n", rank, left);
    }
    else if (!odd_rank && rank != size - 1)
    { // if last rank number is even don't try and receive from the right yet
      MPI_Sendrecv(&data[(work_rows+1) *columns], 1,speed_rows, right, tag, &data[(work_rows+2)*columns], 1, speed_rows, right, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      // printf("rank %d 4\n", rank);

      // printf("round 1 rank= %d send to %d\n", rank, right);
    }
    // printf("round 1 rank=%d cells[0]=%f cells[columns]=%f cells[20]=%f\n",rank,cells[0],cells[columns],cells[20]);

    if (!odd_rank && !(rank == 0 && even_size))
    { // even rank send left
      MPI_Sendrecv(&data[columns], 1, speed_rows, left, tag, data, 1, speed_rows, left, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else if (odd_rank)
    {
      MPI_Sendrecv(&data[(work_rows+1) *columns], 1,speed_rows, right, tag, &data[(work_rows+2)*columns], 1, speed_rows, right, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    // printf("round 2 rank=%d cells[0]=%f cells[columns]=%f cells[20]=%f\n", rank, cells[0], cells[columns], cells[20]);

    if (even_size)
    {
      if (rank == 0)
      {
        MPI_Sendrecv(&data[columns], 1, speed_rows, left, tag, data, 1, speed_rows, left, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      // printf("round 3 rank= %d send to %d\n", rank, left);
      }
      else if (rank == size - 1)
      {
        MPI_Sendrecv(&data[(work_rows+1) *columns], 1,speed_rows, right, tag, &data[(work_rows+2)*columns], 1, speed_rows, right, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      // printf("round 3 rank= %d send to %d\n", rank, right);
      }
    }
  }
  if (rank == 0)
  {
    printf("top halo at rank 0  = %f\n", cells->speed0[0]);
    printf("bottom halo at rank 0  = %f\n", cells->speed0[((work_rows + 2) * columns) );
  }
  else if (rank == 3)
  {
    printf("top halo at rank 3  = %f\n", cells->speed0[0]);
    printf("bottom halo at rank 3 = %f\n", cells->speed0[((work_rows + 2) * columns)]);
  }

  float *global;


/*MPI_Datatype local_cells;
  MPI_Type_vector(9, (work_rows + 1) * columns, (work_rows + 3) * columns, MPI_FLOAT, &local_cells);
  MPI_Type_commit(&local_cells); */

/*   if (rank==0) {
    printf("global size %d\n", (rows + size * 2) * columns * 9);
    global = (float*)_mm_malloc(sizeof(float) * (rows+ size*2) * columns * 9, 64);
  }
  MPI_Gather(data, (work_rows + 3) * columns * 9, MPI_FLOAT, global, (work_rows + 3) * columns * 9, MPI_FLOAT, 0, MPI_COMM_WORLD);

 */
  MPI_Finalize();

/*   if (rank==0) {
    int current_node = 0;
    int work_len = (work[current_node].end_row +1)* columns;
    work_rows = work[current_node].end_row - work[current_node].start_row;
    int node_offset = 0;
    int jj_offset = 0;

    for (size_t jj = 0; jj < rows; jj++)
    {
      for (size_t ii = 0; ii < columns; ii++)
      {
          if (jj > work[current_node].end_row)
          {
            jj_offset = jj;
            node_offset += (work_rows + 3) * columns * 9;
            current_node += 1;
            work_rows = work[current_node].end_row - work[current_node].start_row;
          }
          int index = ii + (jj - jj_offset) *columns; 
          int baseIndex = node_offset + columns + index;
          int speedOffset = (work_rows + 3) * columns;
          printf("speed0[%d]=%f\n", index, global[baseIndex]);
          printf("speed1[%d]=%f\n", index, global[baseIndex + speedOffset*8]);
      }
      
    } */
  
/*     for (size_t n = 0; n < size; n++)
    {
      for (size_t x = 0; x < (work_rows+1)*columns; x++)
      {
        if (rank==0) {
          printf("speed0[%d]=%f\n", x, global[baseIndex + speed]);
          printf("speed1[%d]=%f\n", x, global[node_offset + (work_rows + 3) * columns * 3 + columns + x]);
        }
      }
      node_offset += (work_rows + 3) * columns * 9;
      if (rank==0) {
        //printf("node_offset=%d\n",node_offset);
      }
      current_node += 1;

    }  
    
  }*/
  

  /*   _mm_free(cells->speed0);
    _mm_free(cells->speed1);
    _mm_free(cells->speed2);
    _mm_free(cells->speed3);
    _mm_free(cells->speed4);
    _mm_free(cells->speed5);
    _mm_free(cells->speed6);
    _mm_free(cells->speed7);
    _mm_free(cells->speed8); */
  _mm_free(cells);
  _mm_free(data);
  if (rank==0) {
    //_mm_free(global);
  }
 

  /*   _mm_free(top_data);
    _mm_free(top_halo);
    _mm_free(bottom_data);
    _mm_free(bottom_halo);
   */

  return EXIT_SUCCESS;
}

