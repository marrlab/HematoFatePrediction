function  ss = inprod_basis(basis1, basis2, Lfdobj1, Lfdobj2, rng, ...
                            wtfd, Mmat, EPS, JMAX, JMIN)
%  INPROD_BASIS  Computes matrix of inner products of bases.
%    If both are B-spline bases, both Lfdobj's are
%       numeric, there is no wtfd argument, and the ranges are the same,
%       the inner products are exact, and computed by inprod_bspline.
%    Otherwise, by numerical integration using Romberg integration 
%       with the trapezoidal rule.
%
%  Arguments:
%  BASIS1 and BASIS2 ...  these are basis objects.
%  Lfdobj1 and Lfdobj2 ...  differential operators for inner product for
%                    BASIS1 and BASIS2, respectively
%  RNG  ...  Limits of integration
%  WTFD ...  A functional data object defining a weight
%  MMAT ...  A metric matrix for the vector inner product computed
%            at each step.  The number of rows must equal to
%            NBASIS1 and number columns to NBASIS2.  Defaults to empty.
%  EPS  ...  A convergence criterion, defaults to 1e-4.
%  JMAX ...  Maximum number of Richardson extrapolation iterations.
%            Defaults to 15.
%  JMIN ...  Minimum number of Richardson extrapolation iterations.
%            Defaults to 5.
%
%  Return:
%  A NREP1 by NREP2 matrix SS of inner products for each possible pair
%  of basis functions.

%  last modified 27 October 2008

%  set up default values of arguments

if nargin < 10, JMIN = 5;             end
if nargin <  9, JMAX = 15;            end
if nargin <  8, EPS = 1E-4;           end
if nargin <  7, Mmat = [];            end
if nargin <  6, wtfd    = 0;          end
if nargin <  4, Lfdobj2 = int2Lfd(0); end
if nargin <  3, Lfdobj1 = int2Lfd(0); end
if nargin <  2, basis2  = basis1;     end

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
end
 
%  check BASIS1 and BASIS2

if ~(isa_basis(basis1) && isa_basis(basis2))
    error ('The two first arguments are not basis objects.');
end

%  determine NBASIS1 and NASIS2, and check for common range

nbasis1 = getnbasis(basis1) - length(getdropind(basis1));
nbasis2 = getnbasis(basis2) - length(getdropind(basis2));
type1   = getbasistype(basis1);
type2   = getbasistype(basis2);
range1  = getbasisrange(basis1);
if nargin < 5, rng = range1; end
if rng(1) < range1(1) || rng(2) > range1(2)
    error('Limits of integration are inadmissible.');
end

%  call B-spline version if both bases are the same B-spline basis
%  else proceed with the use of the Romberg integration
  
if strcmp(type1,'bspline')      && strcmp(type2,'bspline')     && ...
    basis1 == basis2                                           && ...
    isinteger(Lfdobj1)          && isinteger(Lfdobj2)          && ...
    nargin < 6                                                 && ...
    isempty(getdropind(basis1)) && isempty(getdropind(basis2)) && ...
    all(rng == range1)
 
    fdobj1 = fd(eye(nbasis1), basis1);
    fdobj2 = fd(eye(nbasis2), basis2);
    ss     = inprod_bspline(fdobj1, fdobj2, ...  
                            getnderiv(Lfdobj1), getnderiv(Lfdobj2));
elseif strcmp(type1,'fourier')  && strcmp(type2,'fourier')     && ...
    getbasispar(basis1) == getbasispar(basis2)                 && ...
    isinteger(Lfdobj1)          && isinteger(Lfdobj2)          && ...
    nargin < 6                                                 && ...
    isempty(getdropind(basis1)) && isempty(getdropind(basis2)) && ...
    all(rng == range1)

    nbasis1 = getnbasis(basis1);
    nbasis2 = getnbasis(basis2);
    ss = zeros(nbasis1, nbasis2);
    if nbasis1 < nbasis2
        ss(1:nbasis1,1:nbasis1) = eval_penalty(basis1, Lfdobj1);
    else
        ss(1:nbasis2,1:nbasis2) = eval_penalty(basis2, Lfdobj2);
    end

else
 
  %  check the metric matrix
  
  if ~isempty(Mmat)
      [m1,m2] = size(Mmat);
      if (nbasis1 ~= m1) || (nbasis2 ~= m2)
          error('MMAT has wrong dimensions.');
      end
  end
  
  %  set up first iteration

    width = rng(2) - rng(1);
    JMAXP = JMAX + 1;
    h     = ones(JMAXP,1);
    h(2)  = 0.25;
    s = reshape(zeros(JMAXP*nbasis1*nbasis2,1),[JMAXP,nbasis1,nbasis2]);
    %  the first iteration uses just the endpoints
    basismat1 = eval_basis(rng, basis1, Lfdobj1);
    basismat2 = eval_basis(rng, basis2, Lfdobj2);
    if ~isnumeric(wtfd)
        wtd = eval_fd(wtfd, rng);
        basismat2 = (wtd * ones(1,nbasis2)) .* basismat2;
    end
    if isempty(Mmat)
        chs = width.*(basismat1' * basismat2)./2;
    else
        chs = width.*(basismat1' * Mmat * basismat2)./2;
    end
    s(1,:,:) = chs;
    tnm = 0.5;

    %  now iterate to convergence

    for iter = 2:JMAX
        tnm = tnm.*2;
        del = width./tnm;
        x   = rng(1)+del/2:del:rng(2);
        basismat1 = eval_basis(x, basis1, Lfdobj1);
        basismat2 = eval_basis(x, basis2, Lfdobj2);
        if ~isnumeric(wtfd)
            wtd = eval_fd(wtfd, x);
            basismat2 = (wtd * ones(1,nbasis2)) .* basismat2;
        end
        if isempty(Mmat)
            chs = width.*(basismat1' * basismat2)./tnm;
        else
            chs = width.*(basismat1' * Mmat * basismat2)./tnm;
        end
        chsold = reshape(s(iter-1,:,:),size(chs));
        if nbasis1 == 1 || nbasis2 == 1
            s(iter,:)   = (chsold + chs)'./2;
        else
            s(iter,:,:) = (chsold + chs)./2;
        end
        if iter >= 5
            ind = (iter-4):iter;
            ya = s(ind,:,:);
            xa = h(ind);
            absxa = abs(xa);
            [absxamin, ns] = min(absxa);
            cs = ya;
            ds = ya;
            y  = squeeze(ya(ns,:,:));
            ns = ns - 1;
            for m = 1:4
                for i = 1:(5-m)
                    ho      = xa(i);
                    hp      = xa(i+m);
                    w       = (cs(i+1,:,:) - ds(i,:,:))./(ho - hp);
                    ds(i,:,:) = hp.*w;
                    cs(i,:,:) = ho.*w;
                end
                if 2*ns < 5-m
                    dy = squeeze(cs(ns+1,:,:));
                else
                    dy = squeeze(ds(ns,:,:));
                    ns = ns - 1;
                end
                y = y + dy;
            end
            ss = reshape(y, nbasis1, nbasis2);
            errval = max(max(abs(dy)));
            ssqval = max(max(abs(ss)));
            if all(ssqval > 0)
                crit = errval./ssqval;
            else
                crit = errval;
            end
            if crit < EPS && iter >= JMIN
                return
            end
        end
        s(iter+1,:,:) = s(iter,:,:);
        h(iter+1)   = 0.25.*h(iter);
    end
    disp(['No convergence after ',num2str(JMAX),' steps in INPROD_BASIS.']);
  
end

