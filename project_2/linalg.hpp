#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <time.h>
#include "jacobi.hpp"
#ifndef linalg_H
#define linalg_H





void DeallocateMatrix(double ** Matrix, int m, int n){
  for(int i=0;i<m;i++)
    delete[] Matrix[i];
  delete[] Matrix;
}


double ** AllocateMatrix(int m, int n){
  double ** Matrix;
  Matrix = new double*[m];
  for(int i=0;i<m;i++){
    Matrix[i] = new double[n];
    for(int j=0;j<m;j++)
      Matrix[i][j] = 0.0;
  }
  return Matrix;
}

double ** MatrixMultiplication(double ** a, double **b, int n){
  double **c = AllocateMatrix(n, n);
  for(int i=0;i < n;i++){
     for (int j=0 ; j < n;j++){
       double sum = 0.0;
       for (int k = 0; k < n; k++) sum += a[i][k]*b[k][j];
       c[i][j] = sum;
     }
  }
  return c;
}