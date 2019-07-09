//
// Created by Daniel Owen on 2019-05-20.
//

#include "helperFuncs.h"
#include <cmath>
#include <stdio.h>

#ifndef C___LAMBERT_BATTIN_H
#define C___LAMBERT_BATTIN_H

void lambert_battin(vector R1, vector R2, double dT, double mu, double dir, vector V[2]);

double battin_xi(double x);
double battin_K(double u);

#endif //C___LAMBERT_BATTIN_H
