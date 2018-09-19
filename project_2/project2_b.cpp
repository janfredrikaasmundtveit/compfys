#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <time.h>
//#include <armadillo>
// use namespace for output and input
using namespace std;
//using namespace arma;


void Jacobi_rotate ( double**, double**, int, int, int);
double offdiag(double**, int*, int*,int);
void DeallocateMatrix(double **, int, int);
double ** AllocateMatrix(int, int);
double ** MatrixMultiplication(double **, double **, int);


int main(int argc, char *argv[])
{
  int n=atoi(argv[1]);
  double **A=AllocateMatrix(n,n);
  double **R=AllocateMatrix(n,n);
  double **C=AllocateMatrix(n,n);
  int maxiter=100;
  double tolerance = 1.0E-10; 
double maxnondiag=1.0;
int iterations = 0;
for (int i = 0; i < n; i++)
{for (int j = 0; j < n; j++)
{
  if (i==j)
  {
    A[i][j]=-2.0;
    R[i][j]=1.0;
    C[i][j]=1.0;
  }
  else if(i==j+1 || i==j-1){
    A[i][j]=1.0;
   
  }
  else {

    A[i][j]=0.0;
     }
}

}
/*for (int i = 0; i < n; i++)
{ cout << endl;
  for (int j = 0; j < n; j++)
{
  
cout << A[i][j] << '|';  
}}*/
double sum, rij;
while ( maxnondiag > tolerance && iterations <= maxiter)
{ for (int i = 0; i < n; i++)
{for (int j = 0; j < n; j++)
{ C[i][j]=R[i][j];
 if (i==j)
  {
    R[i][j]=1.0;
  } 
  else{
    R[i][j]=0.0;
  }
}}
  
   int p, q;
   maxnondiag= offdiag(A, &p, &q, n);
   Jacobi_rotate(A, R, p, q, n);
   iterations++;
   R=MatrixMultiplication(R,C,n);
    cout<< p <<','<< q << endl; 
}

/*
for (int i = 0; i < n; i++)
{ cout << endl;
  for (int j = 0; j < n; j++)
{
  
cout << A[i][j] << '|';  
}}
cout << endl;
for (int i = 0; i < n; i++)
{ cout << endl;
  for (int j = 0; j < n; j++)
{
  
cout << R[i][j] << '|';  
}}

//test R^TBR=A
/*
double **B=AllocateMatrix(n,n);
for (int i = 0; i < n; i++)
{for (int j = 0; j < n; j++)
{C[i][j]=R[j][i];
}
}
B=MatrixMultiplication(A,R,n);
B=MatrixMultiplication(C,B,n);
for (int i = 0; i < n; i++)
{ cout << endl;
  for (int j = 0; j < n; j++)
{
  
cout << B[i][j] << '|';  
}}
DeallocateMatrix(B,n,n);
*/

cout << endl;
cout << iterations;
DeallocateMatrix(A,n,n);
DeallocateMatrix(R,n,n);
DeallocateMatrix(C,n,n);

return 0;
}

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



