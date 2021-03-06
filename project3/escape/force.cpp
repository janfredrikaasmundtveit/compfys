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
#include "totforce.h"


coord force(planet p1, planet p2,double b){
double pi=acos(-1.0);
double c=4*pi*pi;
coord F;
double R=p1.r()+p2.r();
F.x=c*p1.m*p2.m*(p1.q.x-p2.q.x)/(pow(R,b));
F.y=c*p1.m*p2.m*(p1.q.y-p2.q.y)/(pow(R,b));

return F;
}
