//
// Created by Daniel Owen on 2019-06-21.
//

#ifndef C___LAMBERT_BATTIN_MULTI_H
#define C___LAMBERT_BATTIN_MULTI_H

#include "Lambert_Battin.h"

#define pi 4*atan(1)

void lambert_battin_multi(vector R1, vector R2, double dT, double mu, double dir, int N, vector V[2]);

double revSucSub(int N, double m, double l, double x0, double y);
double sucSub(int N, double m, double l, double x0, double y);
double ySolve(double c);

#endif //C___LAMBERT_BATTIN_MULTI_H
