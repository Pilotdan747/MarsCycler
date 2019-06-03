//
// Created by Daniel Owen on 2019-05-20.
//

#include "Lambert_Battin.h"

#define pi 4*atan(1)

double eps = 10-15;

void lambert_battin(vector R1, vector R2, double dT, double mu, double dir, vector V[2]) {
    double r1, r2, kk, theta, c, s, L, T, r0p, l, m, Tp, x0, x, a;
    bool flag = 0;


    r1 = norm(R1);
    r2 = norm(R2);

    kk = cross(R1, R2).z;

    theta = acos(dot(R1, R2)/r1/r2);

    if ((dir == 0 && kk < 0) || (dir > 0 && kk >= 0)) {
        theta = 2 * pi - theta;
        flag = 1;
    }

    c = sqrt(pow(r1, 2) + pow(r2, 2) - 2*r1*r2*cos(theta));
    s = 0.5*(r1 + r2 + c);

    L = sqrt((s-c)/s);
    if (theta > pi)
        L = -L;

    //T = sqrt(8*mu/pow(s, 3))*dT;

    r0p = 0.25*s*pow((1+L), 2);
    l = pow(((1-L)/(1+L)),2);
    //m = pow(T, 2)/pow((1+L), 6);
    m = mu*pow(dT, 2)/(8*pow(r0p, 3));

    //Tp = 4.0/3.0*pow((1-L), 3);
    if (flag)
        Tp = sqrt(2)*pow(s, 3.0/2.0)/(3*sqrt(mu))*(1 - pow((s - c)/s, 3.0/2.0));
    else
        Tp = sqrt(2)*pow(s, 3.0/2.0)/(3*sqrt(mu))*(1 + pow((s - c)/s, 3.0/2.0));

    if (dT <= Tp)
        x0 = 0;
    else
        x0 = l;

    x = x0;


    double z, den, h1, h2, B, u, K, y;

    for (int i = 0; i < 100; i++) {
        z = battin_xi(x);
        if (isnan(z)) {
            int a = 1;
        }

        den = (1 + 2 * x + l) * (4 * x + z * (3 + x));
        h1 = pow((l + x), 2) * (1 + 3 * x + z) / den;
        h2 = m * (x - l + z) / den;
        B = 0.25 * 27 * h2 / pow((1 + h1), 3);
        u = 0.5 * B / (1 + sqrt(1 + B));
        K = battin_K(u); //Check this
        y = (1 + h1) / 3 * (2 + sqrt(1 + B) / (1 + 2 * u * pow(K, 2)));
        x = sqrt(0.25 * pow((1 - l), 2) + m / pow(y, 2)) - 0.5 * (1 + l);
        if (abs(x - x0) < eps)
            break;
        else
            x0 = x;
    }

    a = mu*pow(dT, 2)/16.0/pow(r0p, 2)/x/pow(y, 2);

    double b, amin, tmin, ae, dE, f, g, gdot, ah, bh, dH;

    if (a > 0) {
        b = 2 * asin(sqrt(0.5 * (s - c) / a));
        if (theta > pi) {
            b = -b;
        }

        amin = 0.5 * s;
        tmin = sqrt(pow(amin, 3) / mu) * (pi - b + sin(b));
        ae = 2 * asin(sqrt(0.5 * s / a));

        if (dT > tmin) {
            ae = 2 * pi - ae;
        }

        dE = ae - b;
        f = 1 - a / r1 * (1 - cos(dE));
        g = dT - sqrt(pow(a, 3) / mu) * (dE - sin(dE));
        gdot = 1 - a / r2 * (1 - cos(dE));
    } else {
        ah = 2 * asinh(sqrt(-0.5 * s / a));
        bh = 2 * asinh(sqrt(-0.5 * (s - c) / a));
        dH = ah - bh;
        f = 1 - a / r1 * (1 - cosh(dH));
        g = dT - sqrt(pow(-a, 3) / mu) * (sinh(dH) - dH);
        gdot = 1 - a / r2 * (1 - cosh(dH));
    }

    vector V1, V2;

    V1.x = (R2.x - f*R1.x)/g; V1.y = (R2.y - f*R1.y)/g; V1.z = (R2.z - f*R1.z)/g;
    V2.x = (gdot*R2.x - R1.x)/g; V2.y = (gdot*R2.y - R1.y)/g; V2.z = (gdot*R2.z - R1.z)/g;

    V[0] = V1;
    V[1] = V2;
}


double battin_xi(double x) {
    double tiny, d, n, f0, C0, D0, D, C, Del, f;

    tiny = 1e-30;
    d = sqrt(1+x)+1;
    n = x/(d*d);

    f0 = tiny; C0 = f0; D0 = 0;

    D = 3 + 8*d*D0;
    if (abs(D) < tiny)
        D = tiny;

    C = 3 + 8*d/C0;
    if (abs(C) < tiny)
        C = tiny;

    D = 1.0/D;
    Del = C*D;
    f = f0*Del;
    f0 = f; C0 = C; D0 = D;

    D = 5+n + 1*D0;
    if (abs(D) < tiny)
        D = tiny;

    C = 5+n + 1.0/C0;
    if (abs(C) < tiny)
        C = tiny;

    D = 1.0/D;  Del = C*D;
    f = f0*Del;
    f0 = f; C0 = C; D0 = D;

    D = 1 + 9.0/7.0*n*D0;
    if (abs(D) < tiny)
        D = tiny;

    C = 1 + 9.0/7.0*n/C0;
    if (abs(C) < tiny)
        C = tiny;

    D = 1.0/D;  Del = C*D;
    f = f0*Del;
    f0 = f; C0 = C; D0 = D;


    double c;
    for (int i = 0; i < 100; i++) {
        c = pow((i + 3), 2) / (pow((2 * (i + 3)), 2) - 1);
        D = 1 + c * n * D0;
        if (abs(D) < tiny)
            D = tiny;

        C = 1 + c * n / C0;
        if (abs(C) < tiny)
            C = tiny;

        D = 1 / D;
        Del = C * D;
        f = f0 * Del;
        if (abs(Del - 1) < eps)
            break;
        else

            f0 = f;
            C0 = C;
            D0 = D;
    }

    return f;
}

double battin_K(double u) {
    double tiny, f0, C0, D0, D, C, Del, f;
    tiny = 1e-30;

    f0 = tiny; C0 = f0; D0 = 0;

    D = 1 + 1.0/3.0*D0;
    if (abs(D) < tiny)
        D = tiny;

    C = 1 + 1.0/3.0/C0;
    if (abs(C) < tiny)
        C = tiny;

    D = 1/D;  Del = C*D;
    f = f0*Del;
    f0 = f; C0 = C; D0 = D;

    D = 1 + 4.0/27.0*u*D0;
    if (abs(D) < tiny)
        D = tiny;

    C = 1 + 4.0/27.0*u/C0;
    if (abs(C) < tiny)
        C = tiny;

    D = 1.0/D;  Del = C*D;
    f = f0*Del;
    f0 = f; C0 = C; D0 = D;

    double c1, c2;
    for (int i = 0; i < 100; i++) {
        c1 = 2 * (3 * i + 1) * (6 * i - 1) / 9 / (4 * i - 1) / (4 * i + 1);
        c2 = 2 * (3 * i + 2) * (6 * i + 1) / 9 / (4 * i + 1) / (4 * i + 3);
        D = 1 + c1 * u * D0;
        if (abs(D) < tiny)
            D = tiny;

        C = 1 + c1 * u / C0;
        if (abs(C) < tiny)
            C = tiny;

        D = 1.0 / D;
        Del = C * D;
        f = f0 * Del;

        f0 = f;
        C0 = C;
        D0 = D;
        D = 1 + c2 * u * D0;
        if (abs(D) < tiny)
            D = tiny;

        C = 1 + c2 * u / C0;
        if (abs(C) < tiny)
            C = tiny;

        D = 1.0 / D;
        Del = C * D;
        f = f0 * Del;

        if (abs(Del - 1) < eps)
            break;
        else
            f0 = f;
            C0 = C;
            D0 = D;
    }

    return f;
}
