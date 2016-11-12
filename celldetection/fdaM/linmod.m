function linmodstr = linmod(xfd, yfd, betacell, wtvec)
%  LINMOD  Fits an unrestricted or full functional linear model of the form
%       y(t) = \alpha(t) + \int x(s) \beta(s,t) ds + \epsilon(t),
%  where
%       \beta(s,t) = \phi'(s) B \psi(t)
%  
%  Arguments:
%  XFD      ... a functional data object for the independent variable 
%  YFD      ... a functional data object for the   dependent variable 
%  BETACELL ... a cell array of length 3 containing functional
%               parameter objects for \alpha and for 
%               \beta as a function of s and t, respectively.
%  WTVEC   ... a vector of weights
%  Returns:  a struct object LINMODSTR with fields
%  ALPHA ... a functional parameter object for \alpha
%  BETA  ... a bivariate functional parameter object for \beta
%  YHAT  ... a functional data object for the approximation to y

%  Last modified 1 May 2009

%  check xfd and yfd

if ~isa_fd(xfd)
    error('XFD is not a functional data object.');
end

if ~isa_fd(yfd)
    error('YFD is not a functional data object.');
end

ybasis  = getbasis(yfd);
ynbasis = getnbasis(ybasis);
ranget  = getbasisrange(ybasis);

xbasis  = getbasis(xfd);
ranges  = getbasisrange(xbasis);

nfine = max([201,10*ynbasis+1]);
tfine = linspace(ranget(1),ranget(2),nfine)';

%  get dimensions of data

coefy   = getcoef(yfd);
coefx   = getcoef(xfd);
coefdx  = size(coefx);
coefdy  = size(coefy);
ncurves = coefdx(2);
if coefdy(2) ~= ncurves
    error ('Numbers of observations in first two arguments do not match.');
end

%  set up or check weight vector

if nargin < 4
    wtvec = [];
else
    wtvec = wtcheck(ncurves, wtvec);
end

%  get basis parameter objects

alphafdPar  = betacell{1};
betabifdPar = betacell{2};
alphafdPar  = fdParcheck(alphafdPar);
betabifdPar = bifdParcheck(betabifdPar);

%  get Lfd objects

alphaLfd    = getLfd(alphafdPar);
betaLfdcell = getLfd(betabifdPar);
betasLfd    = betaLfdcell{1};
betatLfd    = betaLfdcell{2};

%  get smoothing parameters

alphalambda = getlambda(alphafdPar);
[betaslambda, betatlambda] = getlambda(betabifdPar);

%  get basis objects

alphabasis = getbasis(getfd(alphafdPar));
alpharange = getbasisrange(alphabasis);
if alpharange(1) ~= ranget(1) || alpharange(2) ~= ranget(2)
    error('Range of alpha coefficient and YFD not compatible.');
end

betabifd   = getfd(betabifdPar);

betasbasis = getsbasis(betabifd);
betasrange = getbasisrange(betasbasis);
if betasrange(1) ~= ranges(1) || betasrange(2) ~= ranges(2)
    error('Ranges of betas coefficient and XFD not compatible.');
end
betatbasis = gettbasis(betabifd);
betatrange = getbasisrange(betatbasis);
if betatrange(1) ~= ranget(1) || betatrange(2) ~= ranget(2)
    error('Range of betat coefficient and YFD not compatible.');
end

%  get numbers of basis functions

alphanbasis = getnbasis(alphabasis);
betasnbasis = getnbasis(betasbasis);
betatnbasis = getnbasis(betatbasis);

%  compute inner products of basis functions

alphattmat = inprod(alphabasis, alphabasis);
betalttmat = inprod(betatbasis, alphabasis);
betassmat  = inprod(betasbasis, betasbasis);
betattmat  = inprod(betatbasis, betatbasis);

%  compute inner products of basis functions and data functions

Fmat = inprod(yfd, alphabasis);
Gmat = inprod(yfd, betatbasis);
Hmat = inprod(xfd, betasbasis);

if isempty(wtvec)
    HHCP = Hmat'*Hmat;
    HGCP = Hmat'*Gmat;
    H1CP = sum(Hmat)';
    F1CP = sum(Fmat)';
else
    HHCP = Hmat'*((wtvec*ones(1,betasnbasis)).*Hmat);
    HGCP = Hmat'*((wtvec*ones(1,betatnbasis)).*Gmat);
    H1CP = Hmat'*wtvec;
    F1CP = Fmat'*wtvec;
end

%  get penalty matrices

if alphalambda > 0
    alphapenmat = eval_penalty(alphabasis, alphaLfd);
else
    alphapenmat = [];
end
if betaslambda > 0
    betaspenmat = eval_penalty(betasbasis, betasLfd);
else
    betaspenmat = [];
end
if betatlambda > 0
    betatpenmat = eval_penalty(betatbasis, betatLfd);
else
    betatpenmat = [];
end

%  set up coefficient matrix and right side for stationary equations

betan = betasnbasis*betatnbasis;
ncoef = alphanbasis + betan;
Cmat  = zeros(ncoef,ncoef);
Dmat  = zeros(ncoef,1);

%  rows for alpha

ind1 = 1:alphanbasis;
ind2 = ind1;
Cmat(ind1,ind2) = ncurves.*alphattmat;
if alphalambda > 0
    Cmat(ind1,ind2) = Cmat(ind1,ind2) + alphalambda.*alphapenmat;
end
ind2 = alphanbasis + (1:betan);
Cmat(ind1,ind2) = kron(H1CP,betalttmat)';

Dmat(ind1) = F1CP;

%  rows for beta

ind1 = alphanbasis + (1:betan);
ind2 = 1:alphanbasis;
Cmat(ind1,ind2) = Cmat(ind2,ind1)';
ind2 = ind1;
Cmat(ind1,ind2) = kron(HHCP,betattmat);
if betaslambda > 0
    Cmat(ind1,ind2) = Cmat(ind1,ind2) + ...
                      betaslambda.*kron(betaspenmat,betattmat);
end
if betatlambda > 0
    Cmat(ind1,ind2) = Cmat(ind1,ind2) + ...
                      betatlambda.*kron(betassmat,betatpenmat);
end

Dmat(ind1) = reshape(HGCP',betan,1);

%  solve the equations

coefvec = symsolve(Cmat, Dmat);

%  set up the coefficient function estimates

%  functional structure for the alpha function

ind1 = 1:alphanbasis;
alphacoef = coefvec(ind1);

alphafdnames = getnames(yfd);
alphafdnames{3} = 'Intercept';
alphafd = fd(alphacoef, alphabasis, alphafdnames);

%  bi-functional structure for the beta function

ind1 = alphanbasis + (1:betan);
betacoef    = reshape(coefvec(ind1),betatnbasis,betasnbasis)';
betafdnames = getnames(xfd);
betafdnames{3} = 'Reg. Coefficient';
betafd = bifd(betacoef, betasbasis, betatbasis, betafdnames);

%  functional data structure for the yhat functions

xbetacoef = (Hmat*betacoef)';
xbetafd   = fd(xbetacoef, betatbasis);
yhatmat   = eval_fd(tfine, alphafd)*ones(1,ncurves) + ...
            eval_fd(tfine, xbetafd);
yhatfd    = smooth_basis(tfine, yhatmat, ybasis); 

linmodstr.alpha = alphafd;
linmodstr.beta  = betafd;
linmodstr.yhat  = yhatfd;


