// tesitng gausian eliminatoin code
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
// use namespace for output and input
using namespace std;
 
 ofstream ofile;
// Functions used
inline double f(double x){return x;}


int main(int argc, char const *argv[])
{
	
	    int n=atoi(argv[1])-2;
	double *b = new double [n+1]; double *c=new double [n+1]; double *a = new double [n+1]; double *u = new double [n+1];
      double *x = new double[n+1]; double *g =new double [n+1];

      u[0] = u[n] = 0.0; // boundary conditions
      b[0] = b[n] = 2; // fixed elements of diagonal
     
      for (int i = 0; i <= n; i++){
		
        g[i] = f(i); 
        b[i]=2; // defining the matix (this is stupid)
        a[i]=-1;
        c[i]=-1;
      }
       	// gaussian eliminaton
	 for (int i = 1; i < n; i++) {
		b[i]=b[i]-(a[i]*c[i]/b[i-1]); // back sub
		g[i]=g[i]-(a[i]*g[i-1])/b[i-1]; 
      } 

       u[n-1] = g[n-1]/b[n-1];
      for (int i = n-2; i > 0; i--){
       u[i] = (g[i]+c[i]*u[i+1])/b[i]; // forward sub
    
       
}

cout << u[0] u[1]  u[2] u[3] endl;
	return 0;
}