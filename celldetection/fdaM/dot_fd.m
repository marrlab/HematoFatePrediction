function  [inprodmat, iter] = ...
                     dot_fd(fdobj1, fdobj2, Lfdobj1, Lfdobj2, ...
                            rng, wtfd, EPS, JMAX, JMIN)
%  INPROD   Computes matrix of inner products of functions.
%  These functions can be:
%  (1) regular functional data objects
%  (2) basis objects, in which case they are converted to functional
%      data objects with the identity matrix as the coefficient matrix.
%  Previous versions of INPROD allowed the functions to be bivariate
%  or cell objects, but this has been removed for being too complex
%  a case to manage.
%
%    If both functions have the same B-spline basis, both Lfdobj's are
%       numeric, there is no wtfd argument, and the ranges are the same,
%       the inner products are exact, and computed by inprod_bspline.
%
%    Otherwise, by numerical integration using Romberg integration 
%       with the trapezoidal rule.
%
%  Arguments:
%  FDOBJ1 and FDOBJ2 ... These may be either functional data or basis 
%     function objects.  In the latter case, a functional data object
%     is created from a basis function object by using the
%     identity matrix as the coefficient matrix.
%     If FDOBJ2 is not provided or is empty, it defaults to a function
%     having a constant basis and coefficient 1 for all replications. 
%     This permits the evaluation of simple integrals of functional data
%     objects.
%  LFDOBJ1 and LFDOBJ2 ...  linear differential operators for inner product 
%     for FD1 and FD2, respectively.  Default to 0.  
%  RNG  ...  Limits of integration, defaults to basis range.
%  WTFD ...  A functional data object defining a weight function, 
%            defaults to [].
%  EPS  ...  A convergence criterion, defaults to 1e-4.
%  JMAX ...  Maximum number of Richardson extrapolation iterations.
%            Defaults to 15.
%  JMIN ...  Minimum number of Richardson extrapolation iterations.
%            Defaults to 5.
%  Return:
%  A NREP1 by NREP2 matrix INPRODMAT of inner products for each possible 
%  pair of functions.

%  last modified 10 July 2012

%  set up default values of arguments

if nargin < 9, JMIN = 5;             end
if nargin < 8, JMAX = 15;            end
if nargin < 7, EPS = 1E-4;           end
if nargin < 6, wtfd = [];            end
if nargin < 4, Lfdobj2 = int2Lfd(0); end
if nargin < 3, Lfdobj1 = int2Lfd(0); end

%  check FDOBJ1 for being a functional data object or
%  a bivariate functional data object or a basis or a cell
%  containing functional data objects.
%  If a basis, convert to a functional data object having
%  the identity matrix as the coefficient matrix.

fdclass = 1;
if  isa_fd(fdobj1) 
    coef1  = getcoef(fdobj1);
else
    fdclass = 0;
end

%  Default FDOBJ2 to a constant function, using a basis that matches 
%  that of FDOBJ1 if possible.

if nargin < 2 || isempty(fdobj2)
    tempfd = fdobj1;
    tempbasis = getbasis(tempfd);
    temptype  = getbasistype(tempbasis);
    temprng   = getbasisrange(tempbasis);
    if strcmp(temptype, 'bspline')
        basis2 = create_bspline_basis(temprng, 1, 1);
    elseif strcmp(temptype, 'fourier')
        basis2 = create_fourier_basis(temprng, 1);
    else
        basis2 = create_constant_basis(temprng);
    end
    fdobj2 = fd(1,basis2);
end

if  isa_fd(fdobj2) 
    coef2  = getcoef(fdobj2);
else
    fdclass = 0;
end

% if either FDOBJ1 or FDOBJ2 fails these tests, error message

if ~fdclass
    error (['The two first arguments are not', ...
            ' functional data objects.']);
end

%  Determine NREP1 and NREP2, and check for common range.
%  This must be done differently dependent on whether the
%  functional data object is univariate or bivariate.

%  FDOBJ1

if isa_fd(fdobj1)
    coefd1 = size(coef1);
    if length(coefd1) > 2
        error('Functional data object must be univariate');
    end
    nrep1     = coefd1(2);
    basisobj1 = getbasis(fdobj1);
    type1     = getbasistype(basisobj1);
    range1    = getbasisrange(basisobj1);
else
    error('FDOBJ1 is not a functional data object.');
end

%  FDOBJ2

if isa_fd(fdobj2)
    coefd2 = size(coef2);
    if length(coefd2) > 2
        error('Functional data object must be univariate');
    end
    nrep2     = coefd2(2);
    basisobj2 = getbasis(fdobj2);
    type2     = getbasistype(basisobj2);
    range2    = getbasisrange(basisobj2);
end

%  check that there are the same number of curves

if nrep1 ~= nrep2
    error('The number of curves is not equal.');
else
    nrep = nrep1;
end

%  check for common range and set default range

if any(range1 ~= range2)
    error('Ranges of FDOBJ1 and FDOBJ2 are incompatible.');
end

if nargin < 5  
    rng = range1; 
end

%  check LFDOBJ1 and LFDOBJ2

Lfdobj1 = int2Lfd(Lfdobj1);
Lfdobj2 = int2Lfd(Lfdobj2);
  
%  check WTFD

if  isa_fd(wtfd)
    coefw = getcoef(wtfd);
    coefd = size(coefw);
    if coefd(2) > 1
        error('Argument WTFD is not a single function');
    end
else
    if ~isempty(wtfd)
        error('WTFD is neither empty nor a FD object.');
    end
end
  
%  set iter

iter = 0;

% check range

if rng(1) < range1(1) || rng(2) > range1(2)
    error('Limits of integration are inadmissible.');
end
  
%  Call B-spline version if ...
%  (1) both functional data objects are univariate
%  (2) both bases are B-splines
%  (3) the two bases are identical
%  (4) both differential operators are integers
%  (5) there is no weight function
%  (6) RNG is equal to the range of the two bases.

%  Else proceed with the use of the Romberg integration.
  
if isa_fd(fdobj1)               && isa_fd(fdobj2)          && ...
   strcmp(type1,'bspline')      && strcmp(type2,'bspline') && ...
   basisobj1 == basisobj2       && ...
   isinteger(Lfdobj1)           && isinteger(Lfdobj2)      && ...
   nargin < 6                   && all(rng == range1)
    deriv1 = getnderiv(Lfdobj1);
    deriv2 = getnderiv(Lfdobj2);
    inprodmat = dot_bspline(fdobj1, fdobj2, deriv1, deriv2);
    iter = 0;
    return;
end

%  ------------------------------------------------------------
%  Now determine the number of subintervals within which the
%  numerical integration takes.  This is important if either
%  basis is a B-spline basis and has multiple knots at a 
%  break point.
%  ------------------------------------------------------------

% The default case, no multiplicities.  

rngvec = rng;  

%  check for any knot multiplicities in either argument

knotmult = [];

%  check first functional object for knot multiplicities

%  univariate case
if strcmp(type1,'bspline')
    % Look for knot multiplicities in first basis
    params1  = getbasispar(basisobj1);
    nparams1 = length(params1);
    for i=2:nparams1
        if params1(i) == params1(i-1)
            knotmult = [knotmult, params1(i)];
        end
    end
end

%  check second functional object for knot multiplicities

%  univariate case
if strcmp(type2,'bspline')
    % Look for knot multiplicities in first basis
    params2  = getbasispar(basisobj2);
    nparams2 = length(params2);
    for i=2:nparams2
        if params2(i) == params2(i-1)
            knotmult = [knotmult, params2(i)];
        end
    end
end

%  Modify RNGVEC defining subinvervals if there are any
%  knot multiplicities.

if isempty(knotmult)
    knotmult = sort(unique(knotmult));
    knotmult = knotmult(knotmult > rng(1) & knotmult < rng(2));
    rngvec = [rng(1), knotmult, rng(2)];
end

inprodmat = zeros(nrep, 1);

%  -----------------------------------------------------------------
%                   loop through sub-intervals
%  -----------------------------------------------------------------

nrng = length(rngvec);
for irng = 2:nrng
    rngi = [rngvec(irng-1),rngvec(irng)];
    %  change range so as to avoid being exactly on
    %  multiple knot values
    if irng > 2
        rngi(1) = rngi(1) + 1e-10;
    end
    if irng < nrng
        rngi(2) = rngi(2) - 1e-10;
    end
    
    %  set up first iteration
    
    iter  = 1;
    width = rngi(2) - rngi(1);
    JMAXP = JMAX + 1;
    h     = ones(JMAXP,1);
    h(2)  = 0.25;
    s = zeros(JMAXP,nrep,1);
    %  The first iteration uses just the endpoints.  Evaluate
    %  the objects at these endpoints.
    fx1 = eval_fd(rngi, fdobj1, Lfdobj1);
    fx2 = eval_fd(rngi, fdobj2, Lfdobj2);
    %  multiply by values of weight function if necessary
    if ~isempty(wtfd)
        wtd = eval_fd(rngi, wtfd);
        fx2 = (wtd * ones(1,nrep)) .* fx2;
    end
    
    s(1,:) = width.*fx1.*fx2./2;
    
    tnm = 0.5;
    
    %  now iterate to convergence
    
    for iter = 2:JMAX
        tnm = tnm.*2;
        del = width./tnm;
        x   = rngi(1)+del/2:del:rngi(2);
        fx1 = eval_fd(x, fdobj1, Lfdobj1);
        fx2 = eval_fd(x, fdobj2, Lfdobj2);
        %  multiply by values of weight function if necessary
        if ~isempty(wtfd)
            wtd = eval_fd(x, wtfd);
            fx2 = (wtd * ones(1,nrep)) .* fx2;
        end
        chstemp = width.*fx1.*fx2./tnm;
        s(iter,:) = (s(iter-1,:) + chs)./2;
        if iter >= 5
            ind = (iter-4):iter;
            ya = s(ind,:);
            xa = h(ind);
            absxa = abs(xa);
            [absxamin, ns] = min(absxa);
            cs = ya;
            ds = ya;
            y  = ya(ns,:);
            ns = ns - 1;
            for m = 1:4
                for i = 1:(5-m)
                    ho      = xa(i);
                    hp      = xa(i+m);
                    w       = (cs(i+1,:) - ds(i,:))./(ho - hp);
                    ds(i,:) = hp.*w;
                    cs(i,:) = ho.*w;
                end
                if 2*ns < 5-m
                    dy = squeeze(cs(ns+1,:));
                else
                    dy = squeeze(ds(ns,:));
                    ns = ns - 1;
                end
                y = y + dy;
            end
            ss     = y;
            errval = max(max(abs(dy)));
            ssqval = max(max(abs(ss)));
            if all(ssqval > 0)
                crit = errval./ssqval;
            else
                crit = errval;
            end
            if crit < EPS && iter >= JMIN
                break
            end
        end
        s(iter+1,:) = s(iter,:);
        h(iter+1)   = 0.25.*h(iter);
        if iter == JMAX
            warning('Wid:converge','Failure to converge.');
        end
    end
    
    inprodmat = inprodmat + ss;
    
end

