#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <time.h>

#include <armadillo>
// use namespace for output and input
using namespace std;
inline double V(double x){return 0.0;
}
using namespace arma;
 ofstream ofile;

void Jacobi_rotate ( mat*, mat*, int, int, int);
double offdiag(mat, int*, int*,int);


int main(int argc, char *argv[])
{   string filename;
  filename = argv[1];
  int n=atoi(argv[2]); //segmentaition fault incomming
  double h=1.0/(n);
  double hh=h*h;
  string fileout = filename;
  mat A(n,n);
  mat R(n,n);
 // mat C(n,n);
  int maxiter=100;
  double tolerance = 1.0E-10; 
double maxnondiag=1.0;
int iterations = 0;
for (int i = 0; i < n; i++)
{for (int j = 0; j < n; j++)
{
  if (i==j)
  {
    A(i,j)=-2.0/(hh)+V(i*h);
    R(i,j)=1.0;
    //C(i,j)=1.0;
  }
  else if(i==j+1 || i==j-1){
    A(i,j)=1.0/hh;
   
  }
  else {

    A(i,j)=0.0;
     }
}

}
/*for (int i = 0; i < n; i++)
{ cout << endl;
  for (int j = 0; j < n; j++)
{
  
cout << A[i][j] << '|';  
}}*/

while ( maxnondiag > tolerance && iterations <= maxiter)
{ /*C=R;
  for (int i = 0; i < n; i++)
{for (int j = 0; j < n; j++)
{ 
 if (i==j)
  {
    R(i,j)=1.0;
  } 
  else{
    R(i,j)=0.0;
  }
}}*/
  
   int p, q;
   maxnondiag= offdiag(A, &p, &q, n);
   Jacobi_rotate(&A, &R, p, q, n);
   iterations++;
   //R=R*C;
    //cout<< p <<','<< q << endl; 
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
*/
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
ofile.open(fileout);
      ofile << setiosflags(ios::showpoint | ios::uppercase);
      //      
      for (int i = 0; i < n;i++){for(int j=0; j<n; j++) {
 
          ofile << setw(15) << setprecision(8) << R(i,j);}
         ofile << setw(15) << setprecision(8) <<A(i,i) << endl;}
//cout << endl;
         ofile.close();
cout << iterations;


return 0;
}



void Jacobi_rotate ( mat *B, mat *S, int k, int l, int n )
{ mat A= *B; mat R=*S;
  double s, c;
  if ( A(k,l) != 0.0 ) {
    double t, tau;
    tau = (A(l,l) - A(k,k))/(2*A(k,l));
    
    if ( tau >= 0 ) {
      t = 1.0/(tau + sqrt(1.0 + tau*tau));
    } else {
      t = -1.0/(-tau +sqrt(1.0 + tau*tau));
    }
    
    c = 1/sqrt(1+t*t);
    s = c*t;
  } else {
    c = 1.0;
    s = 0.0;
  }
  double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
  a_kk = A(k,k);
  a_ll = A(l,l);
  A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
  A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
  A(k,l) = 0.0;  // hard-coding non-diagonal elements by hand
  A(l,k) = 0.0;  // same here
  for ( int i = 0; i < n; i++ ) {
    if ( i != k && i != l ) {
      a_ik = A(i,k);
      a_il = A(i,l);
      A(i,k) = c*a_ik - s*a_il;
      A(k,i) = A(i,k);
      A(i,l) = c*a_il + s*a_ik;
      A(l,i) = A(i,l);
    }
//  And finally the new eigenvectors
    r_ik = R(i,k);
    r_il = R(i,l);

    R(i,k) = c*r_ik - s*r_il;
    R(i,l) = c*r_il + s*r_ik;
  } 
  B=A; S=R;
  return;
}
//  the offdiag function, using Armadillo
double offdiag(mat A, int *p, int *q, int n)
{
   double max;
   for (int i = 0; i < n; ++i)
   {
       for ( int j = i+1; j < n; ++j)
       {
           double aij = fabs(A(i,j));
           if ( aij > max)
           { 
              max = aij;  *p = i; *q = j;
           }
       }
   }
return max;
}




