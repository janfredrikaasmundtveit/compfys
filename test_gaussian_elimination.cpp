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
		
	    int n=atoi(argv[2]);
    string filename;
	double *b = new double [n]; double *c=new double [n]; double *a = new double [n]; double *u = new double [n];
      double *x = new double[n]; double *g =new double [n];

       // boundary conditions
      // fixed elements of diagonal
       filename = argv[1];

// Declare new file name
      string fileout = filename;
      // Convert the power 10^i to a string
      string argument = to_string(1);
      // Final filename as filename-i-
      fileout.append(argument);

      for (int i = 0; i < n; i++){
		
        g[i] = f(i); 
       // x[i]=g[i];
        b[i]=2; // defining the matix (this is stupid)
        a[i]=-1;
        c[i]=-1;
      }
       	// gaussian eliminaton
	 for (int i = 1; i < n; i++) {
		b[i]=b[i]-((a[i]*c[i-1])/b[i-1]); // back sub
		g[i]=g[i]-((a[i]*g[i-1])/b[i-1]); 
		     	

      } 

       u[n-1] = g[n-1]/b[n-1];
      for (int i = n-2; i >= 0; i--){
       u[i] = (g[i]-(c[i]*u[i+1]))/b[i]; // forward sub
    	}




      ofile.open(fileout);
      ofile << setiosflags(ios::showpoint | ios::uppercase);
        //  ofile << "       x:             approx:          exact:       relative error" << endl;
      for (int i = 0; i <= n;i++) {
        ofile << setw(15) << setprecision(8) << u[i];
		//ofile << setw(15) << setprecision(8) << x[i];
	      
      }
      ofile.close();
	
 	delete [] a; delete [] b; delete [] c;  delete [] u; delete [] x; delete [] g;
 	

	return 0;
}