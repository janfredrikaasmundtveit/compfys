#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <time.h>
#include <armadillo>

using namespace std;
//verlet 
ofstream ofile;
int main(int argc, char *argv[])
{
 string filename;



if( argc <= 1 ){
    cout << "Bad Usage: " << argv[0] <<
      " read also output file, number of steps, final time on same line" << endl;
    exit(1);
  }
  else{
    filename=argv[1];
   }
  
int n=atoi(argv[2]);
double tfinal=atof(argv[3]);
double pi=acos(-1.0);
double x=1.0; double y=0.0; 
double vx=0.0; double vy=2*pi;	double c=4*pi*pi;
double t=0.0;
double h=tfinal/n;
double r=sqrt(x*x+y*y);
double hh=h*h;
ofile.open(filename);
 ofile << setiosflags(ios::showpoint | ios::uppercase);
		 ofile << setw(15) << setprecision(8) <<t;
	     ofile << setw(15) << setprecision(8) <<x;
	     ofile << setw(15) << setprecision(8) <<y << endl;
	     clock_t tStartj = clock();
while(t<=tfinal){
double prvx=x; double prvy=y;
x=x+h*vx-c*hh*x/(2*r*r*r);
y=y+h*vy-c*hh*y/(2*r*r*r);
vy=vy-(prvy+y)*c*h/(2*r*r*r);
vx=vx-(prvx+x)*c*h/(2*r*r*r);
t=t+h;
r=sqrt(x*x+y*y);
		//t,x,y
	     ofile << setw(15) << setprecision(8) <<t;
	     ofile << setw(15) << setprecision(8) <<x;
	     ofile << setw(15) << setprecision(8) <<y << endl;


}
cout << (double)(clock() - tStartj)/CLOCKS_PER_SEC << endl;
         ofile.close();


	return 0;
}