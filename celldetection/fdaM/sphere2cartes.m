function [coord, Dcoord] = sphere2cartes(angle, radius)
%  SPHERE2CARTES converts angles and radius as polar coordinates for 
%  a point on the n-spnere in n + 1 dimensional space to 
%  cartesian coordinates in COORD, and their derivative of the coordinates
%  with respect to the angles in DCOORD.
%  ANGLE is a vector of n  angles, and RADIUS is the radius of the sphere,

%  Last modified 28 November 2011 by Jim Ramsay

%  radius defaults to 1, check that it is not negative

if nargin < 2, radius = 1;            end
if radius < 0; error('RADIUS < 0.');  end

% S and C are vectors of sines and cosines of angles

angle = angle(:);
n = length(angle);
S = sin(angle);
C = cos(angle);

%  compute coordinates

cumS  = cumprod([1;S]);  %  cumulative products of sines
coord = zeros(n+1,1);
coord(1:n) = radius.*cumS(1:n).*C;
coord(n+1) = radius.*cumS(n+1);

%  compute coordinate derivatives if required

if nargout == 2
    Dcoord = zeros(n+1,n);
    Dcoord(1,1) = -radius*S(1);
    for j=2:n
        for p=1:j-1
            Djp = 1;
            if p > 1
                Djp = Djp*cumS(p);
            end
            Djp = Djp*C(p);
            if p < j-1
                Djp = Djp*prod(S(p+1:j-1));
            end
            if p < j
                Djp = Djp*C(j);
            else
                Djp = Djp*S(j);
            end
            Dcoord(j,p) = radius*Djp;
        end
        Dcoord(j,j) = -radius*cumS(j+1);
    end
    for p=1:n
        Dp = 1;
        if p > 1
            Dp = Dp*cumS(p);
        end
        Dp = Dp*C(p);
        if p < n
            Dp = Dp*prod(S(p+1:n));
        end
        Dcoord(n+1,p) = radius*Dp;
    end
end

