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

% Set Up
phi = deg2rad(30);
dT1 = 120*24*3600;
dT2 = 23*30*24*3600;
dT3 = 120*24*3600;
dT4 = SynodicT*2 - (dT1 + dT2 + dT3); %12*30*24*3600;

%% Earth to Mars CR2BP
Re = re*[1; 0; 0];
Rm = rm*[cos(phi); sin(phi); 0];

dThetaM = 2*pi*(dT1/Tm);
thetaM2 = phi + dThetaM;

dThetaE = 2*pi*(dT1/Te);
thetaE2 = 0 + dThetaE;

Rm2 = rm*[cos(thetaM2); sin(thetaM2); 0];

%Velocities
[V1, V2] = lambert_battin(Re, Rm2, dT1, muSun, 0);

VinfE1 = V1 - Ve*[1; 0; 0];
VinfM2 = V2 - Vm*[cos(thetaM2); sin(thetaM2); 0];
vinfE1 = norm(VinfE1);
vinfM2 = norm(VinfM2);


%% Mars to Mars

dThetaM2 = 2*pi*(dT2/Tm);
thetaM3 = thetaM2 + dThetaM2;

dThetaE2 = 2*pi*(dT2/Te);
thetaE3 = thetaE2 + dThetaE2;

Rm3 = rm*[cos(thetaM3); sin(thetaM3); 0];

%Velocities
[V3, V4] = lambert_battin(Rm2, Rm3, dT2, muSun, 0);

VinfM3 = V3 - Vm*[cos(thetaM2); sin(thetaM2); 0];
VinfM4 = V4 - Vm*[cos(thetaM3); sin(thetaM3); 0];
vinfM3 = norm(VinfM3);
vinfM4 = norm(VinfM4);


%% Mars to Earth

dThetaM3 = 2*pi*(dT3/Tm);
thetaM4 = thetaM3 + dThetaM3;

dThetaE3 = 2*pi*(dT3/Te);
thetaE4 = thetaE3 + dThetaE3;

Re4 = re*[cos(thetaE4); sin(thetaE4); 0];

[V5, V6] = lambert_battin(Rm3, Re4, dT3, muSun, 0);

VinfM5 = V5 - Vm*[cos(thetaM3); sin(thetaM3); 0];
VinfE6 = V6 - Ve*[cos(thetaE4); sin(thetaE4); 0];
vinfM5 = norm(VinfM5);
vinfE6 = norm(VinfE6);

%% Earth to Earth

dThetaM4 = 2*pi*(dT4/Tm);
thetaM5 = thetaM4 + dThetaM4;

dThetaE4 = 2*pi*(dT4/Te);
thetaE5 = thetaE4 + dThetaE4;

Re5 = re*[cos(thetaE5); sin(thetaE5); 0];

[V7, V8] = lambert_battin(Re4, Re5, dT4, muSun, 0);

VinfE7 = V7 - Ve*[cos(thetaE4); sin(thetaE4); 0];
VinfE8 = V8 - Ve*[cos(thetaE5); sin(thetaE5); 0];
vinfE7 = norm(VinfE7);
vinfE8 = norm(VinfE8);

%% Totals

dT = dT1 + dT2 + dT3 + dT4;
dV = abs(vinfM3 - vinfM2) + abs(vinfM5 - vinfM4) + abs(vinfE7 - vinfE6) + abs(vinfE8 - vinfE1);