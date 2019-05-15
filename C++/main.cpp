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

    double re, rm, muSun, Ve, Vm, Te, Tm, SynodicT, phi, dT1, dT2, dT3, dT4, dThetaM, thetaM, dThetaE, thetaE, dT, dV;
    vector Re1, Rm1, Rm2, Rm3, Re4, Re5, vEarth, vMars;

    re = 1.495979e8;
    rm = 2.279483e8;

    muSun = 1.32712440e11;

    Ve = sqrt(muSun/re);
    Vm = sqrt(muSun/rm);

    Te = 2*pi/sqrt(muSun)*pow(re, 3.0 / 2.0);
    Tm = 2*pi/sqrt(muSun)*pow(rm,3.0/2.0);

    SynodicT = 1/(fabs(1/Te - 1/Tm));

    //Set Up & Initial Conditions
    phi = 30*pi/180;

    dT1 = 120*24*3600;
    dT2 = 24*30*24*3600;
    dT3 = 120*24*3600;
    dT4 = SynodicT*2 - (dT1 + dT2 + dT3);

    thetaE = 0;
    thetaM = phi;

    // Earth to Mars
    Re1.x = re*1; Re1.y = re*0; Re1.z = re*0;
    Rm1.x = rm*cos(phi); Rm1.y = rm*sin(phi); Rm1.z = rm*0;

    dThetaM = 2*pi*(dT1/Tm);
    thetaM += dThetaM;

    dThetaE = 2*pi*(dT1/Te);
    thetaE += dThetaE;

    Rm2.x = rm*cos(thetaM); Rm2.y = rm*sin(thetaM); Rm2.z = 0;

    vector V12[2];
    lambert(Re1, Rm2, dT1, muSun, 1, V12);

    vEarth.x = Ve*1; vEarth.y = Ve*0; vEarth.z = Ve*0;
    vMars.x = Vm*cos(thetaM); vMars.y = Vm*sin(thetaM); vMars.z = Vm*0;

    vector VinfE1 = vinf(V12[0], vEarth);
    vector VinfM2 = vinf(V12[1], vMars);

    // Mars to Mars
    dThetaM = 2*pi*(dT2/Tm);
    thetaM += dThetaM;

    dThetaE = 2*pi*(dT2/Te);
    thetaE += dThetaE;

    Rm3.x = rm*cos(thetaM); Rm3.y = rm*sin(thetaM); Rm3.z = rm*0;

    vector V34[2];
    lambert(Rm2, Rm3, dT2, muSun, 1, V34);

    vector VinfM3 = vinf(V34[0], vMars);

    vMars.x = Vm*cos(thetaM); vMars.y = Vm*sin(thetaM); vMars.z = Vm*0;

    vector VinfM4 = vinf(V34[1], vMars);

    // Mars to Earth
    dThetaM = 2*pi*(dT3/Tm);
    thetaM += dThetaM;

    dThetaE = 2*pi*(dT3/Te);
    thetaE += dThetaE;

    Re4.x = re*cos(thetaE); Re4.y = re*sin(thetaE); Re4.z = re*0;

    vector V56[2];
    lambert(Rm3, Re4, dT3, muSun, 1, V56);

    vEarth.x = Ve*cos(thetaE); vEarth.y = Ve*sin(thetaE); vEarth.z = Ve*0;

    vector VinfM5 = vinf(V56[0], vMars);
    vector VinfE6 = vinf(V56[1], vEarth);

    // Earth to Earth
    dThetaM = 2*pi*(dT4/Tm);
    thetaM += dThetaM;

    dThetaE = 2*pi*(dT4/Te);
    thetaE += dThetaE;

    Re5.x = re*cos(thetaE); Re5.y = re*sin(thetaE); Re5.z = re*0;

    vector V78[2];
    lambert(Re4, Re5, dT4, muSun, 1, V78);

    vector VinfE7 = vinf(V78[0], vEarth);

    vEarth.x = re*cos(thetaE); vEarth.y = re*sin(thetaE); vEarth.z = re*0;

    vector VinfE8 = vinf(V78[1], vEarth);

    // Totals
    dT = dT1 + dT2 + dT3 + dT4;
    dV = fabs(norm(VinfM3) - norm(VinfM2)) + fabs(norm(VinfM5) - norm(VinfM4)) + fabs(norm(VinfE7) - norm(VinfE6))
            + fabs(norm(VinfE8) - norm(VinfE1));

    printf("Total Time in Years: %f\n", dT/3600/24/365);
    printf("Total Delta V in Km/s: %f\n", dV);

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
    printf("A is %f\n", A);


    double Znew, Zold, dZ, C, S, y, F, FPri;
    Zold = 0.1;
    int count = 0;
    do {
        C = 0.5 - Zold/24 + pow(Zold, 2)/720;
        S = 1.0/6.0 - Zold/120 + pow(Zold, 2)/5040 - pow(Zold, 3)/362880;
        y = r1 + r2 + A*(Zold*S - 1)/sqrt(C);

        F = pow(y/C, 3.0/2.0)*S + A*sqrt(y) - sqrt(mu)*dT;
        FPri = pow(y/C, 3.0/2.0)*(1/(2*Zold)*(C - 3*S/(2*C)) + 3*pow(S, 2)/(4*C)) + A/8*(3*S/C*sqrt(y) + A*sqrt(C/y));

        Znew = Zold - F/FPri;
        dZ = fabs(Znew - Zold);
        Zold = Znew;
        count++;
        if (count > 1000) {
            printf("Count too high\n");
        }
    } while (dZ > 0.00001 && count < 1001);

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
