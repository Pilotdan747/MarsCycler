function tisserand(num, varargin)
% FUNCTION tisserand
% Author: Brian Kaplinger, Ph.D.
% Last Changed: December 3, 2018
% 
% Arguments:
%   varargin = [vi_step, start, endd, colorp, linew]
%   num - planet number
%   vi_step - interval in excess velocity in km/s
%         (optional, default = 1 )
%   start - excess velocity at which to start contours
%         (optional, default = 1 km/s)
%   endd - excess velocity at which to end contours
%         (optional, default = 10 km/s)
%   colorp - color to draw contours
%         (optional, default = 'b' (blue))
%   linew - line width to draw contours
%         (optional, default = 0.5 pt)
%   
% Use cases:
% tisserand(num) - contours for planet [num], 0.5 pt blue line, 1 km/s to
%     10 km/s in 1 km/s intervals.
% tisserand(num, vi_step) - contours for planet [num], 0.5 pt blue line, 
%     1 km/s to 10 km/s in [vi_step] km/s intervals.
% tisserand(num, vi_step, start) - contours for planet [num], 0.5 pt blue 
%     line, [start] km/s to 10 km/s in [vi_step] km/s intervals.
% tisserand(num, vi_step, start, endd) - contours for planet [num], 0.5 pt 
%     blue line, [start] km/s to [endd] km/s in [vi_step] km/s intervals.
% tisserand(num, vi_step, start, endd, colorp) - contours for planet [num],
%     0.5 pt [colorp] line, [start] km/s to [endd] km/s in [vi_step] km/s 
%     intervals.
% tisserand(num, vi_step, start, endd, colorp, linew) - contours for 
%     planet [num], [linew] pt [colorp] line, [start] km/s to [endd] km/s 
%     in [vi_step] km/s intervals.
%

nargs = length(varargin);
if (nargs == 0)
    vi_step = 1;
    start = 1;
    endd = 10;
    colorp = 'b';
    linew = 0.5;
elseif (nargs == 1)
    vi_step = varargin{1};
    start = 1;
    endd = 10;
    colorp = 'b';
    linew = 0.5;
elseif (nargs == 2)
    vi_step = varargin{1};
    start = varargin{2};
    endd = 10;
    colorp = 'b';
    linew = 0.5;
elseif (nargs == 3)
    vi_step = varargin{1};
    start = varargin{2};
    endd = varargin{3};
    colorp = 'b';
    linew = 0.5;
elseif (nargs == 4)
    vi_step = varargin{1};
    start = varargin{2};
    endd = varargin{3};
    colorp = varargin{4};
    linew = 0.5;
else
    vi_step = varargin{1};
    start = varargin{2};
    endd = varargin{3};
    colorp = varargin{4};
    linew = varargin{5};
end

vis = start:vi_step:endd;
phis = 0:0.01:pi;
n = length(vis);
m = length(phis);

% plantary distances
Rs = [5.79092e7; 1.082073e8; 1.495979e8; 2.279438e8; 7.783408e8;
    1.426666e9; 2.870658e9; 4.498396e9];
mu = 1.32712e11;
AU = Rs(3);
day = 24*3600;

Rp = Rs(num);
Vp = sqrt(mu/Rp);

PERIAPSIS = zeros(n,m); PERIOD = zeros(n,m); nmax = zeros(n,1);
for i = 1:n
    vi = vis(i);
    nmax(i) = m;
    for j = 1:m
        phi = phis(j);
        vperp = Vp - vi*cos(phi);
        vr = vi*sin(phi);
        h = vperp*Rp;
        v = sqrt(vperp^2 + vr^2);
        ee = v^2/2 - mu/Rp;
        a = -mu/2/ee;
        e = sqrt(1-h^2/mu/a);
        rp = a*(1-e);
        T = 2*pi/sqrt(mu)*a^(3/2);
        if (a > 0)
            PERIAPSIS(i,j) = rp/AU;
            PERIOD(i,j) = T/day;
        else
            if (nmax(i) >= m)
            nmax(i) = j-1;
            % break;
            end
        end
    end
end
    

fig = gca;
hold on
for i = 1:n
    plot(fig,PERIOD(i,1:nmax(i)),PERIAPSIS(i,1:nmax(i)),'Color',colorp,...
        'LineWidth',linew)
end
xlabel('PERIOD (days)'); ylabel('RADIUS OF PERIAPSIS (AU)')
hold off
end