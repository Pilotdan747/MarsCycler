% LAMBERT SOLVER USING BATTIN'S METHOD
% By: Brian Kaplinger, Ph.D.
% Changed: February 4, 2017
function [V1, V2] = lambert_battin(R1, R2, dt, mu, dir)
 % Input Variables:
 %   R1: initial position vector (km)
 %   R2: final position vector (km)
 %   dt: elapsed time (s)
 %   mu: standard gravitational parameter (km^3/s^2)
 %   dir: desired solution direction:
 %       dir = 0: prograde
 %       dir = 1: retrograde
 % Output Variables:
 %   V1:  initial velocity vector (km/s)
 %   V2:  final velocity vector (km/s)

r1 = norm(R1); r2 = norm(R2);
kk = [0, 0, 1]*cross(R1,R2);
th = acos(dot(R1,R2)/r1/r2);
if ((dir == 0 && kk < 0) || (dir > 0 && kk >= 0))
    th = 2*pi - th;
end

c = sqrt(r1^2 + r2^2 - 2*r1*r2*cos(th));
s = 0.5*(r1 + r2 + c);

L = sqrt((s-c)/s);
if (th > pi) 
    L = -L; 
end

T = sqrt(8*mu/s^3)*dt;

r0p = 0.25*s*(1+L)^2;
l = ((1-L)/(1+L))^2;
m = T^2/(1+L)^6;

Tp = 4/3*(1-L)^3;
if (T <= Tp)
    x0 = 0;
else 
    x0 = l;
end
x = x0;

for i = 1:100
    z = battin_xi(x);
    den = (1 + 2*x + l)*(4*x + z*(3+x));
    h1 = (l+x)^2*(1+3*x+z)/den;
    h2 = m*(x - l + z)/den;
    B = 0.25*27*h2/(1+h1)^3;
    u = 0.5*B/(1 + sqrt(1+B));
    K = battin_K(u);
    y = (1+h1)/3*(2 + sqrt(1+B)/(1+2*u*K^2));
    x = sqrt(0.25*(1-l)^2+m/y^2) - 0.5*(1+l);
    if (abs(x-x0)<eps)
        break;
    else
        x0 = x;
    end
end

a = mu*dt^2/16/r0p^2/x/y^2;

if (a > 0)
    b = 2*asin(sqrt(0.5*(s-c)/a));
    if (th > pi)
        b = -b;
    end
    amin = 0.5*s;
    tmin = sqrt(amin^3/mu)*(pi - b + sin(b));
    ae = 2*asin(sqrt(0.5*s/a));
    if (dt > tmin)
        ae = 2*pi - ae;
    end
    dE = ae - b;
    f = 1 - a/r1*(1 - cos(dE));
    g = dt - sqrt(a^3/mu)*(dE - sin(dE));
    gdot = 1 - a/r2*(1 - cos(dE));
else
    ah = 2*asinh(sqrt(-0.5*s/a));
    bh = 2*asinh(sqrt(-0.5*(s-c)/a));
    dH = ah - bh;
    f = 1 - a/r1*(1 - cosh(dH));
    g = dt - sqrt(-a^3/mu)*(sinh(dH) - dH);
    gdot = 1 - a/r2*(1 - cosh(dH));
end

V1 = (R2 - f.*R1)./g;
V2 = (gdot.*R2 - R1)./g;
end


function f = battin_xi(x)
tiny = 1e-30;
d = sqrt(1+x)+1;
n = x./(d.*d);

f0 = tiny; C0 = f0; D0 = 0;
% stage 1
D = 3 + 8*d*D0;
if (abs(D) < tiny); D = tiny; end
C = 3 + 8*d/C0;
if (abs(C) < tiny); C = tiny; end
D = 1/D;  Del = C*D;
f = f0*Del;
f0 = f; C0 = C; D0 = D;
% stage 2
D = 5+n + 1*D0;
if (abs(D) < tiny); D = tiny; end
C = 5+n + 1/C0;
if (abs(C) < tiny); C = tiny; end
D = 1/D;  Del = C*D;
f = f0*Del;
f0 = f; C0 = C; D0 = D;
% stage 3
D = 1 + 9/7*n*D0;
if (abs(D) < tiny); D = tiny; end
C = 1 + 9/7*n/C0;
if (abs(C) < tiny); C = tiny; end
D = 1/D;  Del = C*D;
f = f0*Del;
f0 = f; C0 = C; D0 = D;


for i = 1:100
    c = (i+3)^2/((2*(i+3))^2-1);
    D = 1 + c*n*D0;
    if (abs(D) < tiny); D = tiny; end
    C = 1 + c*n/C0;
    if (abs(C) < tiny); C = tiny; end
    D = 1/D;  Del = C*D;
    f = f0*Del;
    if (abs(Del-1) < eps)
        break;
    else
        f0 = f; C0 = C; D0 = D;
    end
end


end

function f = battin_K(u)
tiny = 1e-30;

f0 = tiny; C0 = f0; D0 = 0;
% stage 1
D = 1 + 1/3*D0;
if (abs(D) < tiny); D = tiny; end
C = 1 + 1/3/C0;
if (abs(C) < tiny); C = tiny; end
D = 1/D;  Del = C*D;
f = f0*Del;
f0 = f; C0 = C; D0 = D;
% stage 2
D = 1 + 4/27*u*D0;
if (abs(D) < tiny); D = tiny; end
C = 1 + 4/27*u/C0;
if (abs(C) < tiny); C = tiny; end
D = 1/D;  Del = C*D;
f = f0*Del;
f0 = f; C0 = C; D0 = D;


for i = 1:100
    c1 = 2*(3*i+1)*(6*i-1)/9/(4*i-1)/(4*i+1);
    c2 = 2*(3*i+2)*(6*i+1)/9/(4*i+1)/(4*i+3);
    D = 1 + c1*u*D0;
    if (abs(D) < tiny); D = tiny; end
    C = 1 + c1*u/C0;
    if (abs(C) < tiny); C = tiny; end
    D = 1/D;  Del = C*D;
    f = f0*Del;
    
    f0 = f; C0 = C; D0 = D;
    D = 1 + c2*u*D0;
    if (abs(D) < tiny); D = tiny; end
    C = 1 + c2*u/C0;
    if (abs(C) < tiny); C = tiny; end
    D = 1/D;  Del = C*D;
    f = f0*Del;
    
    if (abs(Del-1) < eps)
        break;
    else
        f0 = f; C0 = C; D0 = D;
    end
end


end
   

