#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <time.h>
#include "factorial.h"
#include "catch.hpp"
#include <armadillo>


 int Factorial(int number ) {
    return number <= 1 ? number : Factorial(number-1)*number;
}