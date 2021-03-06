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
phi = 2.5;
dT1 = 150*24*3600;
dT2 = 26*30*24*3600;
dT3 = 120*24*3600;
dT4 = SynodicT*2 - (dT1 + dT2 + dT3); %12*30*24*3600;

figure
theta = linspace(0, 2*pi);
plot(re*cos(theta), re*sin(theta), 'b');
hold on
plot(rm*cos(theta), rm*sin(theta), 'r');
hold on

N = 1000;
phis = linspace(0, 2*pi, N);
dT1 = linspace(30*24*3600, 12*30*24*3600, N);

count = 1;

for i = 1:N
    for j = 1:N
        dV(count) = cycle(dT1(j), dT2, dT3, dT4, phis(i));
        count = count + 1;
        if mod(count, 1000) == 0
            count/N^2*100
        end
    end
end

figure
plot(1:count - 1, dV)
grid on

function dV = cycle(dT1, dT2, dT3, dT4, phi)
    re = 1.495979e8;
    rm = 2.279483e8;
    muSun = 1.32712440e11;
    Ve = sqrt(muSun/re);
    Vm = sqrt(muSun/rm);
    Te = 2*pi/sqrt(muSun)*re^(3/2);
    Tm = 2*pi/sqrt(muSun)*rm^(3/2);
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
    
    VinfE1 = V1 - Ve*[0; 1; 0];
    VinfM2 = V2 - Vm*[-sin(thetaM2); cos(thetaM2); 0];
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

    VinfM3 = V3 - Vm*[-sin(thetaM2); cos(thetaM2); 0];
    VinfM4 = V4 - Vm*[-sin(thetaM3); cos(thetaM3); 0];
    vinfM3 = norm(VinfM3);
    vinfM4 = norm(VinfM4);


%% Mars to Earth

    dThetaM3 = 2*pi*(dT3/Tm);
    thetaM4 = thetaM3 + dThetaM3;

    dThetaE3 = 2*pi*(dT3/Te);
    thetaE4 = thetaE3 + dThetaE3;

    Re4 = re*[cos(thetaE4); sin(thetaE4); 0];

    [V5, V6] = lambert_battin(Rm3, Re4, dT3, muSun, 0);

    VinfM5 = V5 - Vm*[-sin(thetaM3); cos(thetaM3); 0];
    VinfE6 = V6 - Ve*[-sin(thetaE4); cos(thetaE4); 0];
    vinfM5 = norm(VinfM5);
    vinfE6 = norm(VinfE6);

%% Earth to Earth
    
    dThetaM4 = 2*pi*(dT4/Tm);
    thetaM5 = thetaM4 + dThetaM4;

    dThetaE4 = 2*pi*(dT4/Te);
    thetaE5 = thetaE4 + dThetaE4;

    Rm5 = rm*[cos(thetaM5); sin(thetaM5); 0];
    Re5 = re*[cos(thetaE5); sin(thetaE5); 0];

    [V7, V8] = lambert_battin(Re4, Re5, dT4, muSun, 0);

    VinfE7 = V7 - Ve*[-sin(thetaE4); cos(thetaE4); 0];
    VinfE8 = V8 - Ve*[-sin(thetaE5); cos(thetaE5); 0];
    vinfE7 = norm(VinfE7);
    vinfE8 = norm(VinfE8);

%% Totals

    dT = dT1 + dT2 + dT3 + dT4;
    dV = abs(vinfM3 - vinfM2) + abs(vinfM5 - vinfM4) + abs(vinfE7 - vinfE6) + abs(vinfE8 - vinfE1);
end