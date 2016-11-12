function basismat = sine(evalarg, nbasis, period, nderiv)
%  SINE  Computes the NDERIV derivative of the sine series basis
%    for NBASIS functions with period PERIOD, these being evaluated
%    at values in vector EVALARG.
%  Returns an N by NBASIS matrix BASISMAT of function values

%  last modified 22 December 2012

evalarg = evalarg(:); %  ensure that EVALARG is a column vector
n       = length(evalarg);
onen    = ones(n,1);
range   = [min(evalarg),max(evalarg)];

%  set default number of basis functions

if nargin < 2
    nbasis = n;
end

%  set default period

if nargin < 3
    period = range(2) - range(1);
end

%  set default order of derivative

if nargin < 4
    nderiv = 0;
end

%  check argument values

if nbasis <= 0,  error('NBASIS not positive');  end
if period <= 0,  error('PERIOD not positive');  end
if nderiv <  0,  error('NDERIV negative');      end

%  set up the basis matrix

basismat = zeros(n,nbasis);

%  set up some constants

omega = 2*pi/period;
fac   = omega .* evalarg;

if nderiv == 0
    %  The sine series itself is required.
    basismat(:,1) = 1/sqrt(2);
    j    = 1:(nbasis-1);
    args = fac * j;
    basismat(:,j+1) = sin(args);
else
    %  A derivative of the sine series is required.
    basismat(:,1) = 0.0;
    if nderiv == floor(nderiv/2)*2
        mval  = nderiv/2;
        ncase = 1;
    else
        mval  = (nderiv-1)/2;
        ncase = 2;
    end
    
    j    = 1:(nbasis-1);
    fac  = onen * (((-1).^mval).*(j.*omega).^nderiv);
    args = (omega.*evalarg) * j;
    if ncase == 1
        basismat(:,j+1) =  fac .*  sin(args);
    else
        basismat(:,j+1) = -fac .* -cos(args);
    end
    
end
basismat = basismat./sqrt(period./2);
