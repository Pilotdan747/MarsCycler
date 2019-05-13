% Project 1
% Author: Daniel Owen
% Created: 4/13/19
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clc
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

re = 1.495979e8;
rm = 2.279483e8;

muSun = 1.32712440e11;

Ve = sqrt(muSun/re);
Vm = sqrt(muSun/rm);

Te = 2*pi/sqrt(muSun)*re^(3/2);
Tm = 2*pi/sqrt(muSun)*rm^(3/2);

SynodicT = 1/(abs(1/Te - 1/Tm));

%% Earth to Mars CR2BP
phi = 0;

%Initial
Re = re*[1; 0; 0];
Rm = rm*[cos(phi); sin(phi); 0];

deltaT1 = 90*3600*24; %3600*24*30*(1:0.01:12);

deltaTheta1 = 2*pi*(deltaT1/Tm);
theta1 = phi + deltaTheta1;

%Final
Rm2 = rm*[cos(theta1); sin(theta1); zeros(1, length(theta1))];

for i=1:length(deltaT1)
    [V1, V2] = lambert_battin(Re, Rm2(:,i), deltaT1(i), muSun, 0);

    VinfE(:, i) = V1 - Ve*[1; 0; 0];
    VinfM(:, i) = V2 - Vm*[cos(theta1(i)); sin(theta1(i)); 0];
    vinfE(i) = norm(VinfE(:,i));
    vinfM(i) = norm(VinfM(:,i));
end

deltaT2 = 3600*24*30*(23:0.01:52);

Rm2 = rm*[cos(theta2); sin(theta2); zeros(1, length(theta2))];



%dT1 + dT2 + dT3 + dT4 = SynodicT
%dV1 + dV2 + dV3 + dV4 = dV