#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <time.h>
#include "planet.h"
#include "verlet.h"
#include "force.h"
#include <armadillo>

planet verletstep(planet p, coord F, coord prevF,double h){
p.q.x=p.q.x+h*p.v.x-h*h*F.x/p.m;
p.q.y=p.q.y+h*p.v.y-h*h*F.y/p.m;
p.v.y=p.v.y-(F.y+prevF.y)*h/(2*p.m);
p.v.x=p.v.x-(F.x+prevF.x)*h/(2*p.m);

return p;
}


