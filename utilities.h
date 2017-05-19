#ifndef SUPPORT_H
#define SUPPORT_H

#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include "omp.h"
#define nrowsM1 3 
#define ncolsM1 392 
#define ncolsM2 3


typedef struct {
  char *array;
  size_t used;
  size_t size;
} Array;

void initArray(Array *a, size_t initialSize) {
  a->array = (int *)malloc(initialSize * sizeof(int));
  a->used = 0;
  a->size = initialSize;
}

void insertArray(Array *a, int element) {
  if (a->used == a->size) {
    a->size *= 2;
    a->array = (int *)realloc(a->array, a->size * sizeof(int));
  }
  a->array[a->used++] = element;
}

void freeArray(Array *a) {
  free(a->array);
  a->array = NULL;
  a->used = a->size = 0;
}

void readData(FILE *file, double data[392][4]){
  char c;
  Array ary;
  int i = 0;
  int j = 0;
  initArray(&ary, 1); 
  do{ 
    c = fgetc(file);
   
    if(c == ','){
      data[i][j] = strtod(ary.array,NULL);
      j++;
      freeArray(&ary);
      initArray(&ary, 1); 
    }   
    else if(c == '\n'){
      data[i][j] = strtod(ary.array,NULL);
      j = 0;
      i++;
      freeArray(&ary);
      initArray(&ary, 1); 
    }   
    else{
      insertArray(&ary, c); 
    }   
  }while(c != EOF);
}

double vector_length(double v[392]){
  int i;
  double sum = 0.0;
  double len;
  for(i = 0; i < 392; i++){
    sum += v[i] * v[i];
  }
  len = sqrt(sum);
  return len;
}

double vector_mult(double v1[392], double v2[392]){
  int i;
  double res = 0.0;
  for(i = 0; i < 392; i++){
    res += v1[i] * v2[i];
  }
  return res;
}

void transposeR(double m[3][3], double result[3][3]){
  int i;
  int j;
  #pragma omp parallel for
  for(i = 0; i < 3; i++){
    for(j = 0; j < 3; j++){
      result[j][i] = m[i][j];
    }   
  }
}
void transposeQ(double m[392][3], double result[3][392]){
  int i;
  int j;
  #pragma omp parallel for
  for(i = 0; i < 392; i++){
    for(j = 0; j < 3; j++){
      result[j][i] = m[i][j];
    }
  }
}

void matrix_mult(double m1[nrowsM1][ncolsM1], double m2[ncolsM1][ncolsM2], double result[nrowsM1][ncolsM2]){
  int i, j, k;
  // multiply matrices
  for (k = 0; k < nrowsM1; k++)
    for (i = 0; i < ncolsM2; i++) {
      result[i][k] = 0.0;
      for (j = 0; j < ncolsM1; j++)
        result[k][i] = result[k][i] + m1[k][j] * m2[j][i];
    }
 
}

MPI_Status status;
void matrix_mult_mpi(double m1[nrowsM1][ncolsM1], double m2[ncolsM1][ncolsM2], double result[nrowsM1][ncolsM2]){
  int numtasks,numworkers,source,dest,id;
  int rows,offset;
  int i,j,k;
  int initialized, finalized;

  MPI_Initialized(&initialized);
  if (!initialized){
    MPI_Init(NULL, NULL);
  }
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  numworkers = numtasks-1;

  // master process
  if (id == 0) {
    rows = nrowsM1/numworkers;
    printf("rows: %i\n", rows);
    offset = 0;
    
    // send data to worker processes
    for (dest=1; dest<=numworkers; dest++)
    {
      MPI_Send(&offset, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
      MPI_Send(&rows, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
      MPI_Send(&m1[offset][0], rows*ncolsM1, MPI_DOUBLE,dest,1, MPI_COMM_WORLD);
      MPI_Send(&m2, ncolsM1*ncolsM2, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
      offset = offset + rows;
    }
    
    // recieve results from the worker bees
    for (i = 1; i <= numworkers; i++)
    {
      source = i;
      MPI_Recv(&offset, 1, MPI_INT, source, 2, MPI_COMM_WORLD, &status);
      MPI_Recv(&rows, 1, MPI_INT, source, 2, MPI_COMM_WORLD, &status);
      MPI_Recv(&result[offset][0], rows*ncolsM2, MPI_DOUBLE, source, 2, MPI_COMM_WORLD, &status);
    }
    
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++)
        printf("%6.2f   ", result[i][j]);
      printf ("\n");
    }
  }

  // worker process
  if (id > 0) {
    source = 0;
    MPI_Recv(&offset, 1, MPI_INT, source, 1, MPI_COMM_WORLD, &status);
    MPI_Recv(&rows, 1, MPI_INT, source, 1, MPI_COMM_WORLD, &status);
    MPI_Recv(&m1, rows*ncolsM1, MPI_DOUBLE, source, 1, MPI_COMM_WORLD, &status);
    MPI_Recv(&m2, ncolsM1*ncolsM2, MPI_DOUBLE, source, 1, MPI_COMM_WORLD, &status);

    // multiply matrices
    for (k = 0; k < ncolsM2; k++){
      for (i = 0; i < rows; i++){
        result[i][k] = 0.0;
        for (j = 0; j < ncolsM1; j++){
          result[i][k] = result[i][k] + m1[i][j] * m2[j][k];
        }
      }
    }


    MPI_Send(&offset, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
    MPI_Send(&rows, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
    MPI_Send(&result, rows*ncolsM2, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
  }

  MPI_Finalized(&finalized);
  if (!finalized){
    MPI_Finalize();
  }
}





#endif /* SUPPORT_H */
