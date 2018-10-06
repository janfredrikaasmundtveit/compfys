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


coord totforce(planet p, planet p1,double b){
	coord F,F1;
	F1=force(p,p1,b);
	//F2=force(p,p2);
	F.x=F1.x; F.y=F1.y; 
	return F;
}