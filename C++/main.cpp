#include <iostream>
#include <cmath>

#define pi 4*atan(1)

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

int main() {

    double re = 1.495979e8;
    double rm = 2.279483e8;

    double muSun = 1.32712440e11;

    double Ve = sqrt(muSun/re);
    double Vm = sqrt(muSun/rm);

    double Te = 2*pi/sqrt(muSun)*pow(re, 3.0 / 2.0);
    double Tm = 2*pi/sqrt(muSun)*pow(rm,3.0/2.0);

    double SynodicT = 1/(abs(1/Te - 1/Tm));

    //Set Up & Initial Conditions
    double phi = 30*pi/180;

    double dT1 = 120*24*3600;
    double dT2 = 23*30*24*3600;
    double dT3 = 120*24*3600;
    double dT4 = SynodicT*2 - (dT1 + dT2 + dT3);

    // Earth to Mars
    vector Re1;
    vector Rm1;
    Re1.x = re*1; Re1.y = re*0; Re1.z = re*0;
    Rm1.x = rm*cos(phi); Rm1.y = rm*sin(phi); Rm1.z = rm*0;

    double dThetaM = 2*pi*(dT1/Tm);
    double thetaM2 = phi + dThetaM;

    double dThetaE = 2*pi*(dT1/Te);
    double thetaE2 = 0 + dThetaE;

    vector Rm2;
    Rm2.x = rm*cos(thetaM2); Rm2.y = rm*sin(thetaM2); Rm2.z = 0;

    vector V12[2];
    lambert(Re1, Rm2, dT1, muSun, 1, V12);

    vector vEarth1;
    vEarth1.x = Ve*1; vEarth1.y = Ve*0; vEarth1.z = Ve*0;

    vector vMars2;
    vMars2.x = Vm*cos(thetaM2); vMars2.y = Vm*sin(thetaM2); vMars2.z = Vm*0;

    vector VinfE1 = vinf(V12[0], vEarth1);
    vector VinfM2 = vinf(V12[1], vMars2);

    printf("Test %f", VinfM2.x);
    return 0;
}

double norm(vector v) {
    return sqrt(pow(v.x, 2) + pow(v.y, 2) + pow(v.z, 2));
}

vector vinf(vector v, vector vPlanet) {
    vector vInf;

    vInf.x = v.x - vPlanet.x;
    vInf.y = v.y - vPlanet.y;
    vInf.z = v.z - vPlanet.z;

    return vInf;
}

vector cross(vector a, vector b) {
    vector c;

    c.x = a.y * b.z - a.z - b.y;
    c.y = a.x * b.z - a.z * b.x;
    c.z = a.x * b.y - a.y * b.x;

    return c;
}

double dot(vector a, vector b) {
    double c = 0;

    c += a.x*b.x;
    c += a.y*b.y;
    c += a.z*b.z;

    return c;
}

void lambert(vector R1, vector R2, double dT, double mu, int k, vector V[2]) {

    double r1 = norm(R1);
    double r2 = norm(R2);

    vector c = cross(R1, R2);

    double dTheta;

    if (c.z > 0) {
        dTheta = acos(dot(R1, R2)/ (r1*r2));
    } else {
        dTheta = 2*pi - acos(dot(R1, R2)/ (r1*r2));
    }

    double A = sin(dTheta)*sqrt(r1*r2/(1 - cos(dTheta)));


    double Znew, Zold = 0.1, dZ, C, S, y, F, FPri;
    do {
        C = 0.5 - Zold/24;
        S = 1.0/6.0 - Zold/120;
        y = r1 + r2 + A*(Zold*S - 1)/sqrt(C);

        F = pow(y/C, 3.0/2.0)*S + A*sqrt(y) - sqrt(mu)*dT;
        FPri = pow(y/C, 3.0/2.0)*(1/(2*Zold)*(C - 3*S/(2*C)) + 3*pow(S, 2)/(4*C)) + A/8*(3*S/C*sqrt(y) + A*sqrt(C/y));

        Znew = Zold - F/FPri;
        dZ = abs(Znew - Zold);
        Zold = Znew;
    } while (dZ > 0.00001);

    y = r1 + r2 + A*(Zold*S - 1)/sqrt(C);

    double f = 1 - y/r1;
    double g = A*sqrt(y/mu);
    double gDot = 1 - y/r2;

    vector v1;
    v1.x = 1/g*(R2.x - f*R1.x); v1.y = 1/g*(R2.y - f*R1.y); v1.z = 1/g*(R2.z - f*R1.z);

    vector v2;
    v2.x = 1/g*(gDot*R2.x - R1.x); v2.y = 1/g*(gDot*R2.y - R1.y); v2.z = 1/g*(gDot*R2.z - R1.z);

    V[0] = v1;
    V[1] = v2;
}
