#include <iostream>
#include <cmath>
#include "helperFuncs.h"
#include "cycler.h"

#define pi 4*atan(1)


int main() {
    double dT1, dT2, dT3, dT4, dT, SynodicT, phi, re, rm, muSun, Te, Tm;

    re = 1.495979e8;
    rm = 2.279483e8;

    muSun = 1.32712440e11;

    Te = 2*pi/sqrt(muSun)*pow(re, 3.0 / 2.0);
    Tm = 2*pi/sqrt(muSun)*pow(rm,3.0/2.0);

    SynodicT = 1/(fabs(1/Te - 1/Tm));


    //Set Up & Initial Conditions
    phi = 30*pi/180;

    dT1 = 120*24*3600;
    dT2 = 25*30*24*3600;
    dT3 = 120*24*3600;
    dT4 = SynodicT*2 - (dT1 + dT2 + dT3);
    dT = dT1 + dT2 + dT3 + dT4;

    printf("dT1 is: %4.2f months\n", dT1/3600/24/30);
    printf("dT2 is: %4.2f months\n", dT2/3600/24/30);
    printf("dT3 is: %4.2f months\n", dT3/3600/24/30);
    printf("dT4 is: %4.2f months\n", dT4/3600/24/30);
    printf("Total Time in Months: %f\n", dT/3600/24/30);

    double dV = cycle(dT1, dT2, dT3, phi);

    printf("Total Delta V in Km/s: %f\n", dV);

    return 0;
}
