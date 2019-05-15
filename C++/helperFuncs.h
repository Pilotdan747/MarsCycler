//
// Created by Daniel Owen on 2019-05-15.
//

#ifndef C___HELPERFUNCS_H
#define C___HELPERFUNCS_H

struct vector {
    double x;
    double y;
    double z;
};

double norm(vector v);

vector cross(vector a, vector b);
double dot(vector a, vector b);

vector vinf(vector v, vector vPlanet);

void lambert(vector R1, vector R2, double dT, double mu, int k, vector V[2]);

#endif //C___HELPERFUNCS_H
