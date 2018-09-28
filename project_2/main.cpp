#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <time.h>
#include "linalg.h"
#include "jacobi.h"
#include "catch.hpp"
#include <armadillo>
// use namespace for output and input
using namespace std;
using namespace arma; 
inline double V(double x){return 0.0;
}
//using namespace arma;
 ofstream ofile;

double offdiag(double **A, int* p, int* q, int n);

int main(int argc, char *argv[])
{   string filename;
  filename = argv[1];
  int n=atoi(argv[2]); //segmentaition fault incomming
  double rmax=1.0;
  double rmin=0.0; 
  double h=(rmax-rmin)/(n+1);
  double hh=h*h;
  string fileout = filename;
  double **A=AllocateMatrix(n,n);
  double **R=AllocateMatrix(n,n);
  double **C=AllocateMatrix(n,n);
  mat B(n,n);
  int maxiter=100000;
  double tolerance = 1.0E-10; 
double maxnondiag=1.0;
int iterations = 0;
for (int i = 0; i < n; i++)
{for (int j = 0; j < n; j++)
{
  if (i==j)
  {
    A[i][j]=-2.0/(hh)+V(i);
    R[i][j]=1.0;
    C[i][j]=1.0;
    B(i,j)=A[i][j];  
  }
  else if(i==j+1 || i==j-1){
    A[i][j]=1.0/hh;
    B(i,j)=A[i][j];
  }
  else {

    A[i][j]=0.0;
    //B(i,j)=A[i][j]; 
  }
}

}
/*for (int i = 0; i < n; i++)
{ cout << endl;
  for (int j = 0; j < n; j++)
{
  
cout << A[i][j] << '|';  
}}*/
clock_t tStartj = clock();
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
   maxnondiag=offdiag(A, &p, &q, n);
   Jacobi_rotate(A, R, p, q, n);
   iterations++;
   R=MatrixMultiplication(R,C,n);
    cout<< p <<','<< q << endl;
     
}
 cout <<'Jacobis method' << (double)(clock() - tStartj)/CLOCKS_PER_SEC << endl;
clock_t tStarta = clock();
vec eigen=eig_sym(B);
int k;
double A0=fabs(A[0][0]);
for(int i=0; i < n; i++){
    int j=0;
    while(fabs(A[i][i]-eigen(j)>10E-8)){
        if(fabs(A[i][i])<A0){
          A0=fabs(A[i][i]);
          k=i;
        }

        if(j==n){
          cout<< 'wrong eigenvaues' << endl;
          return 0;
        }
        j++;
    }


}

  
   cout <<'armadillo' << (double)(clock() - tStarta)/CLOCKS_PER_SEC << endl;


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
*/
//test R^TR=I
/*
double **B=AllocateMatrix(n,n);
for (int i = 0; i < n; i++)
{for (int j = 0; j < n; j++)
{C[i][j]=R[j][i];
}
}
B=MatrixMultiplication(C,R,n);

for (int i = 0; i < n; i++)
{ cout << endl;
  for (int j = 0; j < n; j++)
{
  
cout << B[i][j] << '|';  
}}
DeallocateMatrix(B,n,n);
*/
ofile.open(fileout);
      ofile << setiosflags(ios::showpoint | ios::uppercase);
      //rplace on for groundstate      
      for (int i = 0; i < n;i++){
        ofile << setw(15) << setprecision(8) << rmin+i*h;
        for(int j=0; j<n; j++) {
 
          ofile << setw(15) << setprecision(8) << R[i][j];}
         ofile <<  '\n';}
//cout << endl;
         ofile.close();
cout << iterations;
DeallocateMatrix(A,n,n);
DeallocateMatrix(R,n,n);
DeallocateMatrix(C,n,n);

return 0;
}





//  the offdiag function, using Armadillo
double offdiag(double **A, int* p, int* q, int n)
{ //p=new int; q=new int;
  
   double max=0.0;
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