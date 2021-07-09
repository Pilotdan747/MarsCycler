#include "Lambert_Battin.h"
#include "Lambert_Battin_Multi.h"
#include "helperFuncs.h"

int main() {
    double mu, dT, dT2;
    vector R1, R2;

    R1.x = 5738.395737804391; R1.y = 3313.064323937965; R1.z = 0;
    R2.x = 1773.788346646060; R2.y = 6619.868231556691; R2.z = 0;

    mu = 3.986e5;
    dT = 675.0783288822593;

    vector V[2];
    lambert_battin(R1, R2, dT, mu, 0, V);

    dT2 = dT + 5801.1;

    vector V2[2];
    lambert_battin_multi(R1, R2, dT2, mu, 0, 1, V2);

    printf("V1: [%f, %f, %f]\n", V[0].x, V[0].y, V[0].z);
    printf("V2: [%f, %f, %f]\n", V[1].x, V[1].y, V[1].z);

    printf("Multi V1: [%f, %f, %f]\n", V2[0].x, V2[0].y, V2[0].z);
    printf("Multi V2: [%f, %f, %f]\n", V2[1].x, V2[1].y, V2[1].z);

    return 0;
}