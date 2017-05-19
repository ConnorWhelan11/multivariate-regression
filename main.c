#include "utilities.h"

void findQ(double A[392][3], double Q[392][3]){
  int i;
  double a1[392];
  double a2[392];
  double a3[392];
  double q1[392];
  double q2[392];
  for(i = 0; i < 392; i++ ){
    Q[i][0] = A[i][0];
    Q[i][1] = A[i][1];
    Q[i][2] = A[i][2];

    q1[i] = A[i][0];
    a1[i] = A[i][0];
    a2[i] = A[i][1];
    a3[i] = A[i][2];
  }
  
  double len;
  len = vector_length(q1);
  for(i = 0; i < 392; i++ ){
    q1[i] = q1[i] / len;
    Q[i][0] = A[i][0] / len;
  }
  
  double product1;
  product1 = vector_mult(q1, a2);

  for(i = 0; i < 392; i++){
    Q[i][1] = Q[i][1] - product1 * q1[i];
    q2[i] = Q[i][1];
  }

  len = vector_length(q2); 
  for(i = 0; i < 392; i++ ){
    q2[i] = q1[i] / len;
    Q[i][0] = A[i][0] / len;
  }
  
  double product2 = vector_mult(q2, a3);
  double product3 = vector_mult(q1, a3);

  for(i = 0; i < 392; i++ ){
    Q[i][2] = Q[i][2] - product2 * q2[i] - product3 * q1[i];
  }
}

void findR(double A[392][3], double Q[392][3], double R[3][3]){
  double QT[3][392];
  transposeQ(Q, QT);
  matrix_mult(QT, A, R);
  double RT[3][3];
  transposeR(R, RT);
  int i, j;
  for(i = 0; i < 3; i++)
    for(j=0; j < 3; j++)
      R[i][j] = RT[i][j];
}

void back_substitution(double R[3][3],double Q[392][3], double x[3], double b[392]){
  double QT[3][392];
  transposeQ(Q, QT);
  int i, j;
  double temp[3];
  for(i=0; i < 3; i++){
    double sum = 0.0;
    for(j=0; j < 392; j++)
      sum += QT[i][j] * b[j];
    temp[i] = sum;
  }
  x[2] = temp[2] / R[2][2];
  for(i = 1; i >= 0; --i){
    x[i] = temp[i];
    for(j = i + 1; j < 3; ++j){
      x[i] -= R[i][j] * x[j];
    }
    x[i] = x[i] / R[i][i];
  }


  printf("y = ");
  for(i=0; i < 3; i++){
    printf("(%f)*x%i + ", x[i], i);
  }
  printf("error\n");
}


int main(int argc, char **argv){
  FILE *file = fopen("data.csv", "r");
  double data[392][4];
  readData(file, data);
  fclose(file);
  int i;

  // intialize system to be solved
  double A[392][3];
  double x[3];
  double b[392];
  
  for(i = 0; i < 392; i++){
    A[i][0] = data[i][1];
    A[i][1] = data[i][2];  
    A[i][2] = data[i][3];
    b[i] = data[i][0];
  }

  // initialize QR factorization matrices
  double Q[392][3];
  double R[3][3];

  findQ(A, Q);
  findR(A, Q, R);
  back_substitution(R, Q, x, b);

  
}
