#include <iostream>
#include <cmath>
#include "helperFuncs.h"

#define pi 4*atan(1)


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

    dT1 = 150*24*3600;
    dT2 = 28*30*24*3600;
    dT3 = 150*24*3600;
    dT4 = SynodicT*2 - (dT1 + dT2 + dT3);

    printf("dT1 is: %4.2f months\n", dT1/3600/24/30);
    printf("dT2 is: %4.2f months\n", dT2/3600/24/30);
    printf("dT3 is: %4.2f months\n", dT3/3600/24/30);
    printf("dT4 is: %4.2f months\n", dT4/3600/24/30);

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
    lambert(Re4, Re5, dT4, muSun, 1, V78); // 0 for regrograde

    vector VinfE7 = vinf(V78[0], vEarth);

    vEarth.x = Ve*cos(thetaE); vEarth.y = Ve*sin(thetaE); vEarth.z = Ve*0;

    vector VinfE8 = vinf(V78[1], vEarth);

    double test = fabs(norm(VinfE8) - norm(VinfE1));
    printf("Last dV is: %f\n", test);

    // Totals
    dT = dT1 + dT2 + dT3 + dT4;
    dV = fabs(norm(VinfM3) - norm(VinfM2)) + fabs(norm(VinfM5) - norm(VinfM4)) + fabs(norm(VinfE7) - norm(VinfE6))
            + fabs(norm(VinfE8) - norm(VinfE1));

    printf("Total Time in Years: %f\n", dT/3600/24/365);
    printf("Total Delta V in Km/s: %f\n", dV);

    return 0;
}
