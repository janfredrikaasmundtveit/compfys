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
#include "totforce.h"

using namespace std;

ofstream ofile;

void output(planet);

int main(int argc, char *argv[])
{
double pi=acos(-1.0);
int planetnumber=1;
planet earth, sun;

sun.m=1; sun.q.x=0-0; sun.q.y=0.0; sun.v.x=0.0; sun.v.y=0.0; sun.F.x=0.0; sun.F.y=0.0;
earth.m=3*pow(10,-6); earth.q.x=1.0; earth.q.y=0.0; earth.v.x=0.0; earth.v.y=2*pi; earth.F.x=0.0;earth.F.y=0.0;
//jupiter.m=pow(10,-3); jupiter.q.x=10.0; jupiter.q.y=0.0; jupiter.v.x=0.0; jupiter.v.y=2*pi;jupiter.F.x=0.0; jupiter.F.y=0.0;



string filename;

if( argc <= 3 ){
    cout << "Bad Usage: " << argv[0] <<
      " read also output file, number of steps, final time on same line" << endl;
    exit(1);
  }
  else{
    filename=argv[1];
   }
  
int n=atoi(argv[2]);
double tfinal=atof(argv[3]);
double b=atof(argv[4]);
double t=0.0;
double h=tfinal/n;

coord prevF1;
coord prevF2;
prevF1.x=0.0; prevF1.y=0.0;
prevF2.x=0.0; prevF2.y=0.0;
ofile.open(filename);
ofile << setiosflags(ios::showpoint | ios::uppercase);
		 ofile << setw(15) << setprecision(8) <<t;
		 output(earth);
	   //  output(jupiter);
	     ofile  << endl;
while(t<=tfinal){

earth.F=totforce(earth,sun,b);
//jupiter.F=totforce(jupiter,sun,earth);

earth=verletstep(earth,prevF1,h);
//jupiter=verletstep(jupiter,prevF2,h);
prevF1=earth.F; //prevF2=jupiter.F;
t=t+h;
		 ofile << setw(15) << setprecision(8) <<t;
		 output(earth);
	  //   output(jupiter);
	     ofile  << endl;

}

         ofile.close();
	return 0;
}



void output(planet p){

  	ofile << setw(15) << setprecision(8) <<p.q.x;
	ofile << setw(15) << setprecision(8) <<p.q.y; 

}