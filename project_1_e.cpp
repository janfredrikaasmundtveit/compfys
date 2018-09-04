//project 1
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <time.h>
#include <armadillo>
// use namespace for output and input
using namespace std;
 using namespace arma;
 ofstream ofile;
// Functions used
inline double f(double x){return 100.0*exp(-10.0*x);
}
inline double exact(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}


 int main(int argc, char *argv[])
 {
  clock_t tStart = clock();
 	// creating a file containing number of intergrationpoints.
	int n; 
    string filename;
    int exponent;
    if( argc <= 1 ){
          cout << "Bad Usage: " << argv[0] <<
              " read also file name on same line and max power 10^n" << endl;
          exit(1);
    }
        else{
        filename = argv[1]; // first command line argument after name of program
        exponent = atoi(argv[2]);
}
		for (int i = 1; i <= exponent; i++){
      int n = (int) pow(10.0,i);
      // Declare new file name
      string fileout = filename;
      // Convert the power 10^i to a string
      string argument = to_string(i);
      // Final filename as filename-i-
      fileout.append(argument);
      double h = 1.0/(n);
      double hh = h*h;
      double *u = new double [n+1];  double *x = new double [n+1];
      mat A(n+1,n+1,fill::zeros);
      vec b(n+1);
      vec u(n+1);
     // mat L,U,P;
      A(n,n)=2;
      A(0,0)=2;
      A(n,n-1)=-1;
      A(n-1,n)=-1;
      // do LU here.
      b(0)=hh*f(0.0);
      b(n)=hh*f(1.0);
 	for(int i=1;i<=n-1;i++){
 		A(i,i)=2;
 		A(i,i-1)=-1;
 		A(i-1,i)=-1;
 		b(i)=hh*f(i*h);
 	}
 	//lu(L,U,P,A);


 	// making file (not makefile) 
  	 
      ofile.open(fileout);
      ofile << setiosflags(ios::showpoint | ios::uppercase);
      //      ofile << "       x:             approx:          exact:       relative error" << endl;
      for (int i = 1; i < n;i++) {
	double xval = x[i];
 	 double RelativeError = fabs((exact(xval)-u[i])/exact(xval));
         ofile << setw(15) << setprecision(8) << xval;
         ofile << setw(15) << setprecision(8) << u(i);
         ofile << setw(15) << setprecision(8) << exact(xval);
         ofile << setw(15) << setprecision(8) << log10(RelativeError) << endl;
      }
            ofile << setw(15) << setprecision(8) <<  (double)(clock() - tStart)/CLOCKS_PER_SEC << endl;
      ofile.close();

		//	A.print(); 	
	   delete [] u; delete [] x; 
 	}
 	return 0;
 }
