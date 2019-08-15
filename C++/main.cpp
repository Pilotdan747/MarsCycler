#include <iostream>
#include <cmath>
#include "helperFuncs.h"
#include "cycler.h"
#include <time.h>
#include <chrono>
#include <omp.h>

//Defines pi
#define pi 4*atan(1)


int main() {
    double SynodicT, re, rm, muSun, Te, Tm;

    //Radii of Earth and Mars from sun
    re = 1.495979e8;
    rm = 2.279483e8;

    muSun = 1.32712440e11;

    //Orbital Periods of Earth and Mars
    Te = 2*pi/sqrt(muSun)*pow(re, 3.0/2.0);
    Tm = 2*pi/sqrt(muSun)*pow(rm, 3.0/2.0);

    //Synodic Period
    SynodicT = 1/(fabs(1/Te - 1/Tm));

    //Set Up & Initial Conditions
    int dim1, dim2, dim3, dim4; //Number of points in each dimension
    dim1 = 10;
    dim2 = 10;
    dim3 = 10;
    dim4 = 10;

    //Inital phase angle
    //Phi from 0 to 2*pi or 0 to 360

    //Travel time for each leg
    //dT1 70 to 400
    //dT2 23 to 35
    //dT3 70 to 400
    //dT4 Check to see if this is negative
    //dT = dT1 + dT2 + dT3 + dT4;

    //Set up clock for timing
    using namespace std::chrono;
    auto start = high_resolution_clock::now();

    //Main array that stores delta V values

    double ****dV = (double ****) malloc(dim1 * sizeof(double***));

    for (int i = 0; i < dim2; i++) {
        dV[i] = (double ***) malloc(dim2 * sizeof(double**));
        for (int j = 0; j < dim3; j++) {
            dV[i][j] = (double **) malloc(dim3 * sizeof(double*));
            for (int k = 0; k < dim4; k++) {
                dV[i][j][k] = (double *) malloc(dim4 * sizeof(double));
            }
        }
    }

    /*double phi = (0 + 0 * 3.6) * pi / 180;        // 9 -> 8
    double dT1 = (70 + 0 * 3.3) * 24 * 3600;        // 28 -> 27
    double dT2 = (23 + 27 * 0.12) * 30 * 24 * 3600;  // 97 -> 96
    double dT3 = (70 + 0 * 3.3) * 24 * 3600;      // 9 -> 8
    */
    //cycleMulti(dT1, dT2, dT3, phi);

//Main loop region
//Tests all of phi and delta T 1-3 times
//Stores results in the delta V array

#pragma omp parallel
    {
        //Testing code
        int threadID = omp_get_thread_num();
        if (threadID == 0) {
            printf("Num threads is: %d\n", omp_get_num_threads());
            printf("Num procs is: %d\n", omp_get_num_procs());
        }
#pragma omp for
        for (int i = 0; i < dim1; i++) {
            for (int j = 0; j < dim2; j++) {
                for (int k = 0; k < dim3; k++) {
                    for (int l = 0; l < dim4; l++) {
                        //Set dT1-3 and phi for each iteration
                        double phi = (0 + i * 3.6) * pi / 180;        // 9 -> 8
                        double dT1 = (70 + j * 3.3) * 24 * 3600;        // 28 -> 27
                        double dT2 = (23 + k * 0.12) * 30 * 24 * 3600;  // 97 -> 96
                        double dT3 = (70 + l * 3.3) * 24 * 3600;      // 9 -> 8


                        /* double phi = (23.8 + i * 0.1) * pi / 180;          // 28.8
                        double dT1 = (154.1 + j * 0.1) * 24 * 3600;        // 159.1
                        double dT2 = (33.52 + k * 0.02) * 30 * 24 * 3600;  // 34.52
                        double dT3 = (91.4 + l * 0.1) * 24 * 3600;         // 96.4
                        */

                        /*double phi = (29.9 + i * 0.01) * pi / 180;          // 66 30.4
                        double dT1 = (161.5 + j * 0.01) * 24 * 3600;        // 79 162
                        double dT2 = (34.21 + k * 0.005) * 30 * 24 * 3600;  // 47 34.46
                        double dT3 = (96.4 + l * 0.01) * 24 * 3600;         // 55 96.9
                        */
                        double dT4 = SynodicT*2 - (dT1 + dT2 + dT3);

                        double ans;

                        if (dT4 < 0) {
                            ans = 1000;
                        } else {
                            //Calc delta V
                            ans = cycleMulti(dT1, dT2, dT3, phi);
                        }

                        if (isnan(ans )) {
                            //printf("NAN\n");
                        }

                        //Store answer
                        dV[i][j][k][l] = ans;
                    }
                }
            }
        }
    }

    //Stops lock and gets run time
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    double time = duration.count();
    printf("It took %f seconds to run\n", time/1e6);

    //Outputs main array to a CSV file
    //Puts size of each dimentsion in its own line
    //Prints array into 1 huge line
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

    return 0;
}
