//
// Created by Daniel Owen on 2019-05-17.
//

#include "cycler.h"


double cycle(double dT1, double dT2, double dT3, double phi) {
    double re, rm, muSun, Ve, Vm, Te, Tm, SynodicT, dT4, dThetaM, thetaM, dThetaE, thetaE, dT, dV,
            dV1, dV2, dV3, dV4;
    vector Re1, Rm1, Rm2, Rm3, Re4, Re5, vEarth, vMars;

    re = 1.495979e8;
    rm = 2.279483e8;

    muSun = 1.32712440e11;

    Ve = sqrt(muSun/re);
    Vm = sqrt(muSun/rm);

    Te = 2*pi/sqrt(muSun)*pow(re, 3.0 / 2.0);
    Tm = 2*pi/sqrt(muSun)*pow(rm,3.0/2.0);

    SynodicT = 1/(fabs(1/Te - 1/Tm));

    dT4 = SynodicT*2 - (dT1 + dT2 + dT3);

    thetaE = 0;
    thetaM = phi;

// Earth to Mars
    Re1.x = re * 1;
    Re1.y = re * 0;
    Re1.z = re * 0;
    Rm1.x = rm * cos(phi);
    Rm1.y = rm * sin(phi);
    Rm1.z = rm * 0;

    dThetaM = 2 * pi * (dT1 / Tm);
    thetaM += dThetaM;

    dThetaE = 2 * pi * (dT1 / Te);
    thetaE += dThetaE;

    Rm2.x = rm * cos(thetaM);
    Rm2.y = rm * sin(thetaM);
    Rm2.z = 0;

    vector V12[2];
    lambert_battin(Re1, Rm2, dT1, muSun, 0, V12);

    vEarth.x = Ve * 1;
    vEarth.y = Ve * 0;
    vEarth.z = Ve * 0;
    vMars.x = Vm * cos(thetaM);
    vMars.y = Vm * sin(thetaM);
    vMars.z = Vm * 0;

    vector VinfE1 = vinf(V12[0], vEarth);
    vector VinfM2 = vinf(V12[1], vMars);

// Mars to Mars
    dThetaM = 2 * pi * (dT2 / Tm);
    thetaM += dThetaM;

    dThetaE = 2 * pi * (dT2 / Te);
    thetaE += dThetaE;

    Rm3.x = rm * cos(thetaM);
    Rm3.y = rm * sin(thetaM);
    Rm3.z = rm * 0;

    vector V34[2];
    lambert_battin(Rm2, Rm3, dT2, muSun, 0, V34);

    vector VinfM3 = vinf(V34[0], vMars);

    vMars.x = Vm * cos(thetaM);
    vMars.y = Vm * sin(thetaM);
    vMars.z = Vm * 0;

    vector VinfM4 = vinf(V34[1], vMars);

// Mars to Earth
    dThetaM = 2 * pi * (dT3 / Tm);
    thetaM += dThetaM;

    dThetaE = 2 * pi * (dT3 / Te);
    thetaE += dThetaE;

    Re4.x = re * cos(thetaE);
    Re4.y = re * sin(thetaE);
    Re4.z = re * 0;

    vector V56[2];
    lambert_battin(Rm3, Re4, dT3, muSun, 0, V56);

    vEarth.x = Ve * cos(thetaE);
    vEarth.y = Ve * sin(thetaE);
    vEarth.z = Ve * 0;

    vector VinfM5 = vinf(V56[0], vMars);
    vector VinfE6 = vinf(V56[1], vEarth);

// Earth to Earth
    dThetaM = 2 * pi * (dT4 / Tm);
    thetaM += dThetaM;

    dThetaE = 2 * pi * (dT4 / Te);
    thetaE += dThetaE;

    Re5.x = re * cos(thetaE);
    Re5.y = re * sin(thetaE);
    Re5.z = re * 0;

    vector V78[2];
    lambert_battin(Re4, Re5, dT4, muSun, 0, V78);

    vector VinfE7 = vinf(V78[0], vEarth);

    vEarth.x = Ve * cos(thetaE);
    vEarth.y = Ve * sin(thetaE);
    vEarth.z = Ve * 0;

    vector VinfE8 = vinf(V78[1], vEarth);

// Totals
    dV1 = fabs(norm(VinfM3) - norm(VinfM2));
    dV2 = fabs(norm(VinfM5) - norm(VinfM4));
    dV3 = fabs(norm(VinfE7) - norm(VinfE6));
    dV4 = fabs(norm(VinfE8) - norm(VinfE1));
    dV = dV1 + dV2 + dV3 + dV4;

    return dV;
}