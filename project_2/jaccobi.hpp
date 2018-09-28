
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <time.h>
#include "linalg.hpp"

#ifndef jacobi_H
#define jacobi_H


void Jacobi_rotate ( double **A, double **R, int k, int l, int n ){
  double s, c;
  if ( A[k][l] != 0.0 ) {
    double t, tau;
    tau = (A[l][l] - A[k][k])/(2*A[k][l]);
    
    if ( tau >= 0 ) {
      t = 1.0/(tau + sqrt(1.0 + tau*tau));
    } else {
      t = -1.0/(-tau +sqrt(1.0 + tau*tau));
    }
    
    c = 1/sqrt(1+t*t);
    s = c*t;
  }
   else {
    c = 1.0;
    s = 0.0;
  }
  double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
  a_kk = A[k][k];
  a_ll = A[l][l];
  A[k][k] = c*c*a_kk - 2.0*c*s*A[k][l] + s*s*a_ll;
  A[l][l] = s*s*a_kk + 2.0*c*s*A[k][l] + c*c*a_ll;
  A[k][l] = 0.0;  // hard-coding non-diagonal elements by hand
  A[l][k] = 0.0;  // same here
  for ( int i = 0; i < n; i++ ) {
    if ( i != k && i != l ) {
      a_ik = A[i][k];
      a_il = A[i][l];
      A[i][k] = c*a_ik - s*a_il;
      A[k][i] = A[i][k];
      A[i][l] = c*a_il + s*a_ik;
      A[l][i] = A[i][l];
    }

//  And finally the new eigenvectors
    //r_ik = R[i][k];
    //r_il = R[i][l];
  }
    R[l][k] = s;
    R[k][l] =- s;
    R[k][k]=c;
    R[l][l]=c;
  
  return;
}

//  the offdiag function, using Armadillo
double offdiag(double **A, int* p, int* q, int n)
{ //p=new int; q=new int;
  
   double max;
   for (int i = 0; i < n; ++i)
   {
       for ( int j = i+1; j < n; ++j)
       {
           double aij = fabs(A[i][j]);
           if ( aij > max)
           { 
              max = aij;  *p = i; *q = j;
           }
       }
   }
   
   return max;
}
