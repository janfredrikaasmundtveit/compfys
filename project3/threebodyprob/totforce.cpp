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


coord totforce(planet p, planet p1, planet p2){
	coord F,F1,F2;
	F1=force(p,p1);
	F2=force(p,p2);
	F.x=F1.x+F2.x; F.y=F1.y+F2.y; 
	return F;
}