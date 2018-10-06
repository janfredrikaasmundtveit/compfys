#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <time.h>
#include <armadillo>
#include "planet.h"
#include "verlet.h"
#include "force.h"

using namespace std;

ofstream ofile;


int main(int argc, char *argv[])
{
double pi=acos(-1.0);
planet earth;
planet jupiter;



earth.m=3*pow(10,-6); earth.q.x=1.0; earth.q.y=0.0; earth.v.x=0.0; earth.v.y=2*pi;
jupiter.m=pow(10,-3); jupiter.q.x=10.0; jupiter.q.y=0.0; jupiter.v.x=0.0; jupiter.v.y=2*pi;;

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
double t=0.0;
double h=tfinal/n;

coord F;
F=force(earth, jupiter);
coord prevF=F;
ofile.open(filename);
ofile << setiosflags(ios::showpoint | ios::uppercase);
		 ofile << setw(15) << setprecision(8) <<t;
	     ofile << setw(15) << setprecision(8) <<earth.q.x;
	     ofile << setw(15) << setprecision(8) <<earth.q.y; 
	     ofile << setw(15) << setprecision(8) <<jupiter.q.x;
	     ofile << setw(15) << setprecision(8) <<jupiter.q.y << endl;
while(t<=tfinal){

F=force(earth, jupiter);



earth=verletstep(earth,F,prevF,h);
jupiter=verletstep(jupiter,F,prevF,h);
prevF=F;
t=t+h;
		 ofile << setw(15) << setprecision(8) <<t;
	     ofile << setw(15) << setprecision(8) <<earth.q.x;
	     ofile << setw(15) << setprecision(8) <<earth.q.y; 
	     ofile << setw(15) << setprecision(8) <<jupiter.q.x;
	      ofile << setw(15) << setprecision(8) <<jupiter.q.y << endl;


}

         ofile.close();
	return 0;
}