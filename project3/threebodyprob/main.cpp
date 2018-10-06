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
double t=0.0;
double h=tfinal/n;

double pi=acos(-1.0);
int planetnumber=2;
double a=365;
planet earth, jupiter, sun;

sun.m=1; sun.q.x=0.0; sun.q.y=0.0; sun.v.x=0.0; sun.v.y=0.0; sun.F.x=0.0; sun.F.y=0.0;
earth.m=3E-6; earth.q.x=9.763619062330592E-01; earth.q.y= 2.225327099640603E-01; earth.v.x=-3.988919278934853E-03 *a; earth.v.y=1.674541445029773E-02*a; earth.F.x=0.0;earth.F.y=0.0;
jupiter.m=(1.9/2)*1E-3; jupiter.q.x=-2.712032022234562; jupiter.q.y=-4.631713764473182E+00; jupiter.v.x=6.422390077545354E-03*a; jupiter.v.y=-3.452667922526868E-03*a;jupiter.F.x=0.0; jupiter.F.y=0.0;




coord prevF1=totforce(earth,sun,jupiter);
coord prevF2=totforce(jupiter,sun,earth);

ofile.open(filename);
ofile << setiosflags(ios::showpoint | ios::uppercase);
		 ofile << setw(15) << setprecision(8) <<t;
		 output(earth);
	     output(jupiter);
	     ofile  << endl;
while(t<=tfinal){

earth.F=totforce(earth,sun,jupiter);
jupiter.F=totforce(jupiter,sun,earth);

earth=verletstep(earth,prevF1,h);
jupiter=verletstep(jupiter,prevF2,h);
prevF1=earth.F; prevF2=jupiter.F;
t=t+h;
		 ofile << setw(15) << setprecision(8) <<t;
		 output(earth);
	     output(jupiter);
	     ofile  << endl;

}

         ofile.close();
	return 0;
}



void output(planet p){

  	ofile << setw(15) << setprecision(8) <<p.q.x;
	ofile << setw(15) << setprecision(8) <<p.q.y; 

}