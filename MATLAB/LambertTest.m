% Test Lambert Solver
% Author: Daniel Owen
% Date: 2/10/21


%% Test Orbit 200x1000 km Earth orbit
% Single rev -> Theta 30 & 75 degrees

theta1 = deg2rad(30);
theta2 = deg2rad(75);

mu = 3.986e5;

ra = 6378 + 1000;
rp = 6378 + 200;

a = 0.5*(ra+rp);
e = (ra - rp)/(ra + rp);

h = sqrt(mu*a*(1-e^2));

OE1 = [h, e, 0, 0, 0, theta1];
OE2 = [h, e, 0, 0, 0, theta2];

[R1, V1] = OE2SV(mu, OE1)
[R2, V2] = OE2SV(mu, OE2)

E1 = 2*atan(sqrt((1-e)/(1+e))*tan(theta1/2));
E2 = 2*atan(sqrt((1-e)/(1+e))*tan(theta2/2));

Me1 = E1 - e*sin(E1);
Me2 = E2 - e*sin(E2);

T = 2*pi/sqrt(mu)*a^(3/2);

t = T*(Me2 - Me1)/(2*pi)

[Vlam1, Vlam2] = lambert_battin(R1, R2, t, mu, 0)

[VlamMult, VlamMult2] = lambertb(R1, V1, R2, 's', 'd', 0, t) 