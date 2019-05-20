% STATE VECTOR TO ORBITAL ELEMENTS
% By: Brian Kaplinger, Ph.D.
% Changed: September 28, 2018
function OE = SV2OE(mu, R, V)
  % Input Variables:
  %   mu = standard gravitational parameter (km^3/s^2)
  %   R(1:3) = position components  (km)
  %   V(1:3) = velocity components  (km/s)
  % Output Variables:
  %   OE: Vector containing orbital elements (6 elements)
  %   OE(1) = h: specific angular momentum  (km^2/s)
  %   OE(2) = e: eccentricity
  %   OE(3) = i: inclination (radians)
  %   OE(4) = Omega: right ascension of ascending node (radians)
  %   OE(5) = omega: argument of perigee (radians)
  %   OE(6) = theta: true anomaly (radians)
  % Additional Optional Output:
  %   OE(7) = a: energy parameter - semimajor axis for ellipse (km)
  %   OE(8) = p: semiparameter (km)
  
  % Set a small value.  Absolute values lower than this number are
  % considered equal to zero.  In traditional computing theory, this 
  % number is referred to using the greek letter epsilon.  A machine 
  % epsilon is the hardware equivalent of this concept (~1e-14 for double 
  % precision) such that (1 + x) = 1 for all abs(x) < eps.  Matlab will
  % understand "eps" without me setting it as the known machine epsilon
  % (which is a little small for our purposes). Try it!
  eps = 1e-9;
  
  % get magnitudes of state vectors
  r = sqrt(R(1)^2 + R(2)^2 + R(3)^2);
  v = sqrt(V(1)^2 + V(2)^2 + V(3)^2);
  
  % calculate radial velocity vr = dot(R,V)/r
  vr = (R(1)*V(1) + R(2)*V(2) + R(3)*V(3))/r;
  
  % angular momentum vector
  H = cross(R,V);
  
  % get magnitude of angular momentum
  h = sqrt(H(1)^2 + H(2)^2 + H(3)^2);
  
  % calculate eccentricity vector
  E = ((v^2-mu/r).*R - (r*vr).*V)./mu;
  
  % eccentricity is magnitude of this vector
  e = sqrt(E(1)^2 + E(2)^2 + E(3)^2);
  
  % calculate semimajor axis/energy parameter (OPTIONAL)
  if (e > eps)
      a = h^2/mu/(1-e^2);
  else
      a = Inf;
  end
  
  % calculate semiparameter (OPTIONAL)
  p = h^2/mu;
  
  % get ascending node vector 
  Z = [0; 0; 1];
  N = cross(Z,H);
  
  % magnitude of node line vector
  n = sqrt(N(1)^2 + N(2)^2 + N(3)^2);
  
  % calculate inclination 
  i = acos(H(3)/h);
  
  % calculate right ascension of ascending node
  % this is only valid if we are not equatorial (H = Z or n = 0)
  if (n > eps)
      Omega = atan2(N(2),N(1));
      Omega = mod(Omega,2*pi);
  else
      Omega = 0;
  end
 
  % calculate argument of perigee, requires dot(N,E)
  % this is only valid for nonzero inclination (like the right ascension
  % calculation).  For Equatorial orbits, use X axis.  We must also not 
  % have a circular orbit (where would perigee be for a circular orbit?).  
  % So e cannot be zero and still have a well defined argument of perigee.
  if (n < eps)
      N = [1; 0; 0];
      n = 1;
  end
  if (e > eps)
      s = cross(H./h,N./n);
      wy = E(1)*s(1) + E(2)*s(2) + E(3)*s(3);
      wx = (E(1)*N(1) + E(2)*N(2) + E(3)*N(3))/n;
      omega = atan2(wy,wx);
      omega = mod(omega,2*pi);
  else
      omega = 0;
  end
  
  % calculate true anomaly, requires dot(R,E)
  % also if we have a circular orbit handle this case.
  if (e < eps)
      E = N;
      e = n;
  end
  q = cross(H./h,E./e);
  ty = R(1)*q(1) + R(2)*q(2) + R(3)*q(3);
  tx = (R(1)*E(1) + R(2)*E(2) + R(3)*E(3))/e;
  theta = atan2(ty,tx);
  theta = mod(theta,2*pi);
  
  % set output vector
  OE = [h; e; i; Omega; omega; theta; a; p];
end
  