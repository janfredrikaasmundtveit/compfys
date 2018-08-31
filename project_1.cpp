//project 1
//project 1
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
// use namespace for output and input
using namespace std;
 
 ofstream ofile;
// Functions used
inline double f(double x){return 100.0*exp(-10.0*x);
}
inline double exact(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}


 int main(int argc, char *argv[])
 {
 	int n=atoi(argv[1]); 
 	double h=1.0/(n+1.0); //stepsize
	
 	double *b = new double [n+1]; //diagonal elements
 	double *a = new double [n]; //lower off-diagonal
 	double *c = new double [n]; //upper off-diagonal
 	double *x = new double [n+1];
 	double *u = new double [n+1]; // u_n(x)

 	u[0]=0;
 	u[n+1]=0; //boundary conditions

 	//forward sub
 	for(int i=0;i<=n;i++){

 }
 	

 	delete [] a; delete [] b; delete [] c;  delete [] u; delete [] x;
 	return 0;
 }
