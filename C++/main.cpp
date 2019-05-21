#include <iostream>
#include <cmath>
#include "helperFuncs.h"
#include "cycler.h"
#include <time.h>
#include <chrono>


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

    dT1 = 150*24*3600;
    dT2 = 28*30*24*3600;
    dT3 = 100*24*3600;
    dT4 = SynodicT*2 - (dT1 + dT2 + dT3);
    dT = dT1 + dT2 + dT3 + dT4;

    printf("dT1 is: %4.2f months\n", dT1/3600/24/30);
    printf("dT2 is: %4.2f months\n", dT2/3600/24/30);
    printf("dT3 is: %4.2f months\n", dT3/3600/24/30);
    printf("dT4 is: %4.2f months\n", dT4/3600/24/30);
    printf("Total Time in Months: %f\n", dT/3600/24/30);

    using namespace std::chrono;
    auto start = high_resolution_clock::now();

    int dim1, dim2, dim3, dim4;
    dim1 = 20;
    dim2 = 20;
    dim3 = 20;
    dim4 = 20;

    double dV[dim1][dim2][dim3][dim4];

#pragma omp parallel for
    {
        for (int i = 0; i < dim1; i++) {
            for (int j = 0; j < dim2; j++) {
                for (int k = 0; k < dim3; k++) {
                    for (int l = 0; l < dim4; l++) {
                        dT1 = (150 + j * 10) * 24 * 3600;
                        dT2 = (23 + k) * 30 * 24 * 3600;
                        dT3 = (100 + l * 10) * 24 * 3600;
                        phi = i * 10 * pi / 180;
                        dV[i][j][k][l] = cycle(dT1, dT2, dT3, phi);
                    }
                }
            }
        }
    }

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    double time = duration.count();

    printf("It took %f seconds to run\n", time/1e6);

    double dV1 = cycle(dT1, dT2, dT3, 15.5*pi/180);
    printf("dV is %f\n", dV1);

    FILE *outfile = fopen("Output.csv", "w");

    fprintf(outfile, "%d\n", dim1);
    fprintf(outfile, "%d\n", dim2);
    fprintf(outfile, "%d\n", dim3);
    fprintf(outfile, "%d\n", dim4);

    for (int i = 0; i < dim1; i++) {
        for (int j = 0; j < dim2; j++) {
            for (int k = 0; k < dim3; k++) {
                for (int l = 0; l < dim4; l++) {
                    fprintf(outfile, "%f, ", dV[i][j][k][l]);
                }
            }
        }
    }

    fclose(outfile);

    //printf("Total Delta V in Km/s: %f\n", dV);

    return 0;
}
