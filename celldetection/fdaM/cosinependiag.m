function pendiag = cosinependiag(basisobj, nderiv)
%  FOURIERPENDIAG Computes fourier penalty in diagonal case

%  Last modified 29 September 2012

nbasis  = getnbasis(basisobj);
range   = getbasisrange(basisobj);
width   = range(2) - range(1);
params  = getbasispar(basisobj);
period  = params(1);
if width/period ~= round(width/period)
    error('RANGE is not a multiple of PERIOD.');
end
omega   = 2*pi/period;
halfper = period/2;
twonde  = 2*nderiv;
pendiag = zeros(1,nbasis);
if nderiv == 0
    pendiag(1) = period/2;
else
    pendiag(1) = 0;
end
j   = 1:nbasis-1;
fac = halfper*(j*omega).^twonde;
pendiag(j+1) = fac;
pendiag = 2*pendiag./period;
pendiag = (width/period).*pendiag;
