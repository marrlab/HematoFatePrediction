function phimat = IRT3PL(evalarg, nbasis, shift, nderiv)
%  IRT3PL Values of monomials, or their derivatives, combined with
%  log(exp(x) + shift).
%  The powers of EVALARG are the NBASIS nonnegative integers in EXPONENTS.
%  The default is 1, meaning EVALARG itself.
%  Arguments are as follows:
%  EVALARG   ... array of values at which the polynomials are to
%                evaluated
%  SHIFT     ... Shift parameter.   Defaults to 1.
%  NDERIV    ... Order of derivative to be returned.  Max. 2.
%  Return is:
%  A matrix with length(EVALARG) rows and NBASIS columns containing
%    the values of the monomials or their derivatives

%  last modified 19 June 2012

% set default arguments

if nargin < 4,  nderiv = 0;  end
if nargin < 3,  shift  = 1;  end
if nargin < 2,  nbasis = 3;  end

if nderiv > 2
    error('NDERIV exceeds 2.');
end

evalarg = evalarg(:);
n = length(evalarg);

exponents = 0:(nbasis-2);

phimat = zeros(n,nbasis);

if nderiv == 0
    %  use the recursion formula to compute monomnomial values
    for ibasis=1:nbasis-1
        phimat(:,ibasis) = evalarg.^exponents(ibasis); 
    end
else
    for ibasis=1:nbasis-1
        degree = exponents(ibasis);
        if nderiv <= degree 
            fac = degree;
            for ideriv=2:nderiv
                fac = fac*(degree-ideriv+1);
            end
            phimat(:,ibasis) = fac.*evalarg.^(degree-nderiv);
        end
    end
end

%  add log(exp(x) + shift)

evec = exp(evalarg);

switch nderiv
    case 0
        phimat(:,nbasis) = log(evec + shift);
    case 1
        phimat(:,nbasis) = evec./(evec+shift);
    case 2
        phimat(:,nbasis) = evec./(evec+shift).^2;
    otherwise
        error('NDERIV exceeds 2.');
end


