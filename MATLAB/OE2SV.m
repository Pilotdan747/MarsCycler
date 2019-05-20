% ORBITAL ELEMENTS TO STATE VECTOR
% By: Brian Kaplinger, Ph.D.
% Changed: September 28, 2018
function [R, V] = OE2SV(mu, OE)
  % Input Variables:
  %   mu = standard gravitational parameter (km^3/s^2)
  %
  %   OE: Vector containing orbital elements (6 elements)
  %   OE(1) = h: specific angular momentum  (km^2/s)
  %   OE(2) = e: eccentricity
  %   OE(3) = i: inclination (radians)
  %   OE(4) = Omega: right ascension of ascending node (radians)
  %   OE(5) = omega: argument of perigee (radians)
  %   OE(6) = theta: true anomaly (radians)
  % Output Variables:
  %   R(1:3) = position components  (km)
  %   V(1:3) = velocity components  (km/s)
  
  % unroll OE vector for ease of use
  h = OE(1); e = OE(2); i = OE(3); Omega = OE(4); omega = OE(5);
  theta = OE(6);

  % calculate the position vector components in perifocal frame
  r = (h^2/mu)/(1 + e*cos(theta));
  r_p = r.*[cos(theta); sin(theta); 0];
  
  % calculate the velocity vector components in perifocal frame
  v_p = (mu/h).*[-sin(theta); e + cos(theta); 0];
  
  % form rotation matrix of angle Omega about 3 axis
  R1 = [cos(Omega), sin(Omega), 0;
       -sin(Omega), cos(Omega), 0;
                 0,          0, 1];
             
  % form rotation matrix of angle i about 1 axis
  R2 = [1,       0,      0;
        0,  cos(i), sin(i);
        0, -sin(i), cos(i)];
    
  % form rotation matrix of angle omega about 3 axis
  R3 = [cos(omega), sin(omega), 0;
       -sin(omega), cos(omega), 0;
                 0,          0, 1];
  
  % coordinate conversion from geocentric equatorial to perifocal is the
  % matrix multiplication of these three operations (note order)
  C_GP = R3*R2*R1;
  
  % coordinate conversion from perifocal back to geocentric equatorial
  % is the inverse of this matrix, which is also the transpose
  C_PG = C_GP';
  
  % coordinates of position and velocity vectors in geocentric equatiorial
  % frame is this matrix times the array of coordinates in perifocal
  R = C_PG*r_p;
  V = C_PG*v_p;
end
             
             
             