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
    double dT1, dT2, dT3, dT4, dT, SynodicT, phi, re, rm, muSun, Te, Tm;

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
    dim1 = 360;
    dim2 = 370/5;
    dim3 = 12*30/5;
    dim4 = 370/5;

    //Inital phase angle
    phi = 30*pi/180; //0 to 2*pi or 0 to 360

    //Travel time for each leg
    dT1 = 150*24*3600; //70 to 400
    dT2 = 28*30*24*3600; //23 to 35
    dT3 = 100*24*3600; //70 to 400
    dT4 = SynodicT*2 - (dT1 + dT2 + dT3); //Check to see if this is negative
    dT = dT1 + dT2 + dT3 + dT4;

    //Set up clock for timing
    using namespace std::chrono;
    auto start = high_resolution_clock::now();

    //Main array that stores delta V values
    double *dV = (double *) malloc(dim1*dim2*dim3*dim4*sizeof(double));

    //Testing code
    for (int i = 0; i < dim1*dim2*dim3*dim4; i++) {
        dV[i] = rand()*1000;
    }

    //Testing code
    if (dV == NULL) {
        printf("NULL");
    }

    //Testing code
    long test3 = dim1*dim2*dim3*dim4*sizeof(double)/800000;
    printf("Size is: %lu mb\n", test3);


//Main loop region
//Tests all of phi and delta T 1-3 times
//Stores results in the delta V array
#pragma omp parallel
    {
        //Testing code
        int threadID = omp_get_thread_num();
        if (threadID == 0) {
            printf("Num threads is: %d\n", omp_get_num_threads());
        }
#pragma omp for
        for (int i = 0; i < dim1; i++) {
            int count = 0;
            for (int j = 0; j < dim2; j++) {
                for (int k = 0; k < dim3; k++) {
                    for (int l = 0; l < dim4; l++) {
                        //Set dT1-3 and phi for each iteration
                        dT1 = (150 + j * 10) * 24 * 3600;
                        dT2 = (23 + k) * 30 * 24 * 3600;
                        dT3 = (100 + l * 10) * 24 * 3600;
                        phi = (20 + i * 5) * pi / 180;

                        //Calc delta V
                        double ans = cycle(dT1, dT2, dT3, phi);

                        //Testing code
                        /*bool test = isnan(ans);
                        if (isnan(ans)) {
                            printf("dV is: %f Bool showed: %d\n", ans, test);
                        }*/

                        //Store answer
                        dV[count + i*dim2*dim3*dim4] = ans; //Needs work on indexing
                        count++;
                    }
                }
            }
        }
    }

    //Stops lock and gets run time
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    double time = duration.count();
    printf("It took %f seconds to run\n", time/100000);

    //Outputs main array to a CSV file
    //Puts size of each dimentsion in its own line
    //Prints array into 1 huge line
    FILE *outfile = fopen("Output.csv", "w");

    fprintf(outfile, "%d\n", dim1);
    fprintf(outfile, "%d\n", dim2);
    fprintf(outfile, "%d\n", dim3);
    fprintf(outfile, "%d\n", dim4);

    for (int i = 0; i < 1; i++) {
        for (int j = 0; j < 1; j++) {
            for (int k = 0; k < 1; k++) {
                for (int l = 0; l < 1; l++) {
                    //fprintf(outfile, "%f, ", dV[i][j][k][l]);
                }
            }
        }
    }

    fclose(outfile);

    return 0;
}
