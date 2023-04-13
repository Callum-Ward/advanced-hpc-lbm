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

void getHalo(float *buffer, t_speed *cells, int speed_offset, int rank)  {
  for (int i =0;i<10;i++) {
    //printf("rank %d i=%d cells= %f\n", rank,i,cells->speed0[speed_offset*10+i]);
    buffer[i] = cells->speed0[speed_offset*10 + i];
    //printf("rank %d i=%d buffer=%f\n",rank,i,buffer[i]);
    buffer[10+i] = cells->speed1[speed_offset*10 + i];
    buffer[20+i] = cells->speed2[speed_offset*10 + i];
    buffer[30+i] = cells->speed3[speed_offset*10 + i];
    buffer[40+i] = cells->speed4[speed_offset*10 + i];
    buffer[50+i] = cells->speed5[speed_offset*10 + i];
    buffer[60+i] = cells->speed6[speed_offset*10 + i];
    buffer[70+i] = cells->speed7[speed_offset*10 + i];
    buffer[80+i] = cells->speed8[speed_offset*10 + i];
  }
}

void getHalos(float *top_buffer, float *bottom_buffer, t_speed *cells, int rank, int last_row) {
    getHalo(top_buffer,cells,1,rank);
    getHalo(bottom_buffer, cells, last_row,rank);
}

void setHalo(float *buffer, t_speed *cells, int speed_offset) {
    for (int i = 0; i < 10; i++)
    {
      cells->speed0[speed_offset*10 + i] = buffer[i];
      cells->speed1[speed_offset*10 + i] = buffer[10 + i];
      cells->speed2[speed_offset*10 + i] = buffer[20 + i];
      cells->speed3[speed_offset*10 + i] = buffer[30 + i];
      cells->speed4[speed_offset*10 + i] = buffer[40 + i];
      cells->speed5[speed_offset*10 + i] = buffer[50 + i];
      cells->speed6[speed_offset*10 + i] = buffer[60 + i];
      cells->speed7[speed_offset*10 + i] = buffer[70 + i];
      cells->speed8[speed_offset*10 + i] = buffer[80 + i];
    }
}

void setHalos(float *top_buffer, float *bottom_buffer, t_speed *cells, int rank, int last_row)
{
    setHalo(top_buffer, cells, 10);
    setHalo(bottom_buffer, cells, last_row);
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
  get_rank_sizes(size, 8, work);
  // float *cells = (float *)_mm_malloc(sizeof(float) *10* (work[rank].end_row - work[rank].start_row+ 3 ) , 64);

  int work_rows = work[rank].end_row - work[rank].start_row;

  t_speed *cells = (t_speed *)_mm_malloc(sizeof(t_speed), 64);
  cells->speed0 = (float *)_mm_malloc(sizeof(float) * (work_rows+3) * 10, 64);
  cells->speed1 = (float *)_mm_malloc(sizeof(float) * (work_rows+3) * 10, 64);
  cells->speed2 = (float *)_mm_malloc(sizeof(float) * (work_rows+3) * 10, 64);
  cells->speed3 = (float *)_mm_malloc(sizeof(float) * (work_rows+3) * 10, 64);
  cells->speed4 = (float *)_mm_malloc(sizeof(float) * (work_rows+3) * 10, 64);
  cells->speed5 = (float *)_mm_malloc(sizeof(float) * (work_rows+3) * 10, 64);
  cells->speed6 = (float *)_mm_malloc(sizeof(float) * (work_rows+3) * 10, 64);
  cells->speed7 = (float *)_mm_malloc(sizeof(float) * (work_rows+3) * 10, 64);
  cells->speed8 = (float *)_mm_malloc(sizeof(float) * (work_rows+3) * 10, 64);

  t_speed *top_halo = (t_speed *)_mm_malloc(sizeof(t_speed), 64);
  top_halo->speed0 = (float *)_mm_malloc(sizeof(float)* 10, 64);
  top_halo->speed1 = (float *)_mm_malloc(sizeof(float)* 10, 64);
  top_halo->speed2 = (float *)_mm_malloc(sizeof(float)* 10, 64);
  top_halo->speed3 = (float *)_mm_malloc(sizeof(float)* 10, 64);
  top_halo->speed4 = (float *)_mm_malloc(sizeof(float)* 10, 64);
  top_halo->speed5 = (float *)_mm_malloc(sizeof(float)* 10, 64);
  top_halo->speed6 = (float *)_mm_malloc(sizeof(float)* 10, 64);
  top_halo->speed7 = (float *)_mm_malloc(sizeof(float)* 10, 64);
  top_halo->speed8 = (float *)_mm_malloc(sizeof(float)* 10, 64);

  t_speed *bottom_halo = (t_speed *)_mm_malloc(sizeof(t_speed), 64);
  bottom_halo->speed0 = (float *)_mm_malloc(sizeof(float)* 10, 64);
  bottom_halo->speed1 = (float *)_mm_malloc(sizeof(float)* 10, 64);
  bottom_halo->speed2 = (float *)_mm_malloc(sizeof(float)* 10, 64);
  bottom_halo->speed3 = (float *)_mm_malloc(sizeof(float)* 10, 64);
  bottom_halo->speed4 = (float *)_mm_malloc(sizeof(float)* 10, 64);
  bottom_halo->speed5 = (float *)_mm_malloc(sizeof(float)* 10, 64);
  bottom_halo->speed6 = (float *)_mm_malloc(sizeof(float)* 10, 64);
  bottom_halo->speed7 = (float *)_mm_malloc(sizeof(float)* 10, 64);
  bottom_halo->speed8 = (float *)_mm_malloc(sizeof(float)* 10, 64);

  for (int i = 0; i < 10*(work_rows + 3); i++)
  {
    cells->speed0[i] = (float)i;
    cells->speed1[i] = (float)i;
    cells->speed2[i] = (float)i;
    cells->speed3[i] = (float)i;
    cells->speed4[i] = (float)i;
    cells->speed5[i] = (float)i;
    cells->speed6[i] = (float)i;
    cells->speed7[i] = (float)i;
    cells->speed8[i] = (float)i;
  }

  t_speed *cell_trimmed = (t_speed *)_mm_malloc(sizeof(t_speed), 64);
  cell_trimmed->speed0 = &(cells->speed0[10]);
  cell_trimmed->speed1 = &(cells->speed1[10]);
  cell_trimmed->speed2 = &(cells->speed2[10]);
  cell_trimmed->speed3 = &(cells->speed3[10]);
  cell_trimmed->speed4 = &(cells->speed4[10]);
  cell_trimmed->speed5 = &(cells->speed5[10]);
  cell_trimmed->speed6 = &(cells->speed6[10]);
  cell_trimmed->speed7 = &(cells->speed7[10]);
  cell_trimmed->speed8 = &(cells->speed8[10]);


  MPI_Datatype mpi_cells_type;
  const int nitems = 9;
  int blocklengths[9];
  MPI_Datatype types[9];
  MPI_Aint offsets[9];
  for (int i = 0; i < 9; i++)
  {
      blocklengths[i] = 10;
      types[i] = MPI_FLOAT;
  }
  offsets[0] = offsetof(t_speed, speed0);
  offsets[1] = offsetof(t_speed, speed1);
  offsets[2] = offsetof(t_speed, speed2);
  offsets[3] = offsetof(t_speed, speed3);
  offsets[4] = offsetof(t_speed, speed4);
  offsets[5] = offsetof(t_speed, speed5);
  offsets[6] = offsetof(t_speed, speed6);
  offsets[7] = offsetof(t_speed, speed7);
  offsets[8] = offsetof(t_speed, speed8);
  MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_cells_type);
  MPI_Type_commit(&mpi_cells_type);

  printf("created datatype\n");

  if (rank == 0)
  {
    MPI_Sendrecv(cell_trimmed, 1, mpi_cells_type, 1, tag, &bottom_halo, 1, mpi_cells_type, 1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("bottom_halo->speed0[0]=%f\n", bottom_halo->speed0[0]);
  }
  else if (rank == 1)
  {
    MPI_Sendrecv(cell_trimmed, 1, mpi_cells_type, 0, tag, &top_halo, 1, mpi_cells_type, 0, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("top_halo->speed0[0]=%f\n", top_halo->speed0[0]);
  }

  //printf("rank %d cells initialised\n", rank);


  //getHalos(top_Cells,bottom_Cells,cells,rank,work_rows+1);
  //printf("rank %d top cells= %f\n", rank, top_Cells[0]);
  //printf("rank %d bottom cells= %f\n", rank, bottom_Cells[0]);

 /*  for (int i = 0; i <5; i++){

    if (odd_rank) //odds send left 
    {
      MPI_Sendrecv(top_Cells, 90, MPI_FLOAT, left, tag, top_Halo, 90, MPI_FLOAT, left, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      // printf("rank %d 3\n", rank);
      //printf("round 1 rank= %d send to %d\n", rank, left);
    } else if (!odd_rank && rank!=size-1) { //if last rank number is even don't try and receive from the right yet
      MPI_Sendrecv(bottom_Cells, 90,  MPI_FLOAT, right, tag, bottom_Halo, 90,  MPI_FLOAT, right, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      //printf("rank %d 4\n", rank);

      //printf("round 1 rank= %d send to %d\n", rank, right);
    } 
    //printf("round 1 rank=%d cells[0]=%f cells[10]=%f cells[20]=%f\n",rank,cells[0],cells[10],cells[20]);

    if (!odd_rank && !(rank==0 && even_size)) { //even rank send left
      MPI_Sendrecv(top_Cells, 90,  MPI_FLOAT, left, tag, top_Halo, 90,  MPI_FLOAT, left, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      //printf("rank %d 4\n", rank);

      //printf("round 2 rank= %d send to %d\n", rank, left);
    } else if (odd_rank) {
      MPI_Sendrecv(bottom_Cells, 90, MPI_FLOAT, right, tag, bottom_Halo, 90, MPI_FLOAT, right, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      //printf("rank %d 4\n", rank);
      //printf("round 2 rank= %d send to %d\n", rank, right);
    }
    //printf("round 2 rank=%d cells[0]=%f cells[10]=%f cells[20]=%f\n", rank, cells[0], cells[10], cells[20]);

    if (even_size) {
      if (rank==0) {
        MPI_Sendrecv(top_Cells, 90, MPI_FLOAT, left, tag, top_Halo, 90, MPI_FLOAT, left, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //printf("round 3 rank= %d send to %d\n", rank, left);
      } else if (rank==size-1) {
        MPI_Sendrecv(bottom_Cells, 90, MPI_FLOAT, right, tag, bottom_Halo, 90, MPI_FLOAT, right, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //printf("round 3 rank= %d send to %d\n", rank, right);
      }
    } 
  }
 */

  //MPI_Reduce(&cells[10], cells, 10*(work_rows + 1), MPI_FLOAT, MPI_SUM, 0,MPI_COMM_WORLD);
  
  /*   if (rank==0) {
    float *global = (float *)_mm_malloc(sizeof(float) * 80, 64);
    MPI_Gather(&cells[10], 10 * (work_rows + 1), MPI_FLOAT, global, 10 * (work_rows + 1), MPI_FLOAT, 0, MPI_COMM_WORLD);
    for (int i=0;i< 80;i++) {
      printf("total val global[%d] %f\n",i, global[i]);
    }
    _mm_free(global);
  } else {
    MPI_Gather(&cells[10], 9 *10 * (work_rows + 1), MPI_FLOAT, NULL, 10 * (work_rows + 1), MPI_FLOAT, 0, MPI_COMM_WORLD);
  } */

  //printf("rank %d top halo= %f\n", rank, top_Halo[0]);
  //printf("rank %d bottom halo= %f\n", rank, bottom_Halo[0]);


  MPI_Finalize();

 /*  for (int i = 0; i < 10 * (work_rows + 3); i++)
  {
    cells->speed0[i] = (float)i;
    cells->speed1[i] = (float)i;
    cells->speed2[i] = (float)i;
    cells->speed3[i] = (float)i;
    cells->speed4[i] = (float)i;
    cells->speed5[i] = (float)i;
    cells->speed6[i] = (float)i;
    cells->speed7[i] = (float)i;
    cells->speed8[i] = (float)i;
    if (i<9 || i>= 10*(work_rows+2)) {
        cells->speed0[i] = 999.f;
        cells->speed1[i] = 999.f;
        cells->speed2[i] = 999.f;
        cells->speed3[i] = 999.f;
        cells->speed4[i] = 999.f;
        cells->speed5[i] = 999.f;
        cells->speed6[i] = 999.f;
        cells->speed7[i] = 999.f;
        cells->speed8[i] = 999.f;
      }else {
        cells->speed0[i] = (float)rank;
        cells->speed1[i] = (float)rank;
        cells->speed2[i] = (float)rank;
        cells->speed3[i] = (float)rank;
        cells->speed4[i] = (float)rank;
        cells->speed5[i] = (float)rank;
        cells->speed6[i] = (float)rank;
        cells->speed7[i] = (float)rank;
        cells->speed8[i] = (float)rank;
      } 

  }

  t_speed *cell_trimmed = (t_speed *)_mm_malloc(sizeof(t_speed), 64);
  cell_trimmed->speed0 = &(cells->speed0[9]);
  cell_trimmed->speed1 = &(cells->speed1[9]);
  cell_trimmed->speed2 = &(cells->speed2[9]);
  cell_trimmed->speed3 = &(cells->speed3[9]);
  cell_trimmed->speed4 = &(cells->speed4[9]);
  cell_trimmed->speed5 = &(cells->speed5[9]);
  cell_trimmed->speed6 = &(cells->speed6[9]);
  cell_trimmed->speed7 = &(cells->speed7[9]);
  cell_trimmed->speed8 = &(cells->speed8[9]);
  if (rank == 1) {
    printf("rank %d cells_trimmed->speed0[0] = %f\n", rank, cell_trimmed->speed0[0]);
    printf("rank %d cells->speed0[0] = %f\n", rank, cells->speed0[9]);

    printf("&cells->speed0[9]=%p\n", &cells->speed0[9]);
    printf("&cells->speed0[9]=%p\n", &cell_trimmed->speed0[0]);

    for (int i = 0; i < 10 * (work_rows + 3); i++)
    {
      if (i==9) {
        cells->speed0[i] = 5.1;
      } else {
        cells->speed0[i] =(float)i+1;
      } 
      //cells->speed0[i] = (float)i+1;
      cells->speed1[i] =(float)i+1;
      cells->speed2[i] =(float)i+1;
      cells->speed3[i] =(float)i+1;
      cells->speed4[i] =(float)i+1;
      cells->speed5[i] =(float)i+1;
      cells->speed6[i] =(float)i+1;
      cells->speed7[i] =(float)i+1;
      cells->speed8[i] =(float)i+1;
    }

    cells->speed0[9] = 6.1;

    printf("&cells->speed0[9]=%p\n", &cells->speed0[9]);
    printf("&cells->speed0[9]=%p\n", &cell_trimmed->speed0[0]);

    printf("rank %d cells_trimmed->speed0[0] = %f\n", rank, cell_trimmed->speed0[0]);
    printf("rank %d cells->speed0[9] = %f\n", rank, cells->speed0[9]);

    
  }
 */


  _mm_free(cells->speed0);
  _mm_free(cells->speed1);
  _mm_free(cells->speed2);
  _mm_free(cells->speed3);
  _mm_free(cells->speed4);
  _mm_free(cells->speed5);
  _mm_free(cells->speed6);
  _mm_free(cells->speed7);
  _mm_free(cells->speed8);
  _mm_free(cells);
  _mm_free(cell_trimmed);

  _mm_free(top_halo->speed0);
  _mm_free(top_halo->speed1);
  _mm_free(top_halo->speed2);
  _mm_free(top_halo->speed3);
  _mm_free(top_halo->speed4);
  _mm_free(top_halo->speed5);
  _mm_free(top_halo->speed6);
  _mm_free(top_halo->speed7);
  _mm_free(top_halo->speed8);

  _mm_free(bottom_halo->speed0);
  _mm_free(bottom_halo->speed1);
  _mm_free(bottom_halo->speed2);
  _mm_free(bottom_halo->speed3);
  _mm_free(bottom_halo->speed4);
  _mm_free(bottom_halo->speed5);
  _mm_free(bottom_halo->speed6);
  _mm_free(bottom_halo->speed7);
  _mm_free(bottom_halo->speed8);

  _mm_free(top_halo);
  _mm_free(bottom_halo);

  return EXIT_SUCCESS;
}

