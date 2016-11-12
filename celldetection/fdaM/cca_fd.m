function ccastr = ...
              cca_fd(fdobjx, fdobjy, ncan, ccafdParx, ccafdPary, centerfns)
%  CCA_FD   Functional canonical correlation analysis with regularization.
%
%  Arguments:
%  FDOBJx    ... Functional data object for the first  set of functions.
%  FDOBJy    ... Functional data object for the second set of functions.
%                If FDOByx and FDOBJy correspond to two variables, i and j,  
%                of a single multivariate functional data object, say 
%                MVFDOBJ, then use as arguments 
%                MVFDOBJ(:,i) and MVFDOBJ(:,j), respectively.
%  NCAN      ... Number of pairs of canonical variates to be estimated. 
%                Default 2
%  CCAFDPARx ... A functional parameter object for the canonical weight
%                functions for the first  set of functions.
%                The default is fdPar(FDOBJx).
%  CCAFDPARy ... A functional parameter object for the canonical weight
%                functions for the second set of functions.
%                The default is fdPar(FDOBJy).
%  CENTERFNS ... Center the functions before analyzing. Default 1.
%  Returns:  
%  CCASTR   ... A struct object with fields:
%  CCAWTFDx ... A functional data object for the canonical weight
%                functions for the first  set of functions.
%  CCAWTFDy ... A functional data object for the canonical weight
%                functions for the second set of functions.
%  CCAVARx  ... Canonical variate scores for first  set of functions.
%  CCAVARy  ... Canonical variate scores for second set of functions.
%  CORRS    ... The corresponding set of canonical correlations.

%  Last modified on:  24 December 2012

%  check fd objects

if ~strcmp(class(fdobjx), 'fd')
    error('FDOBJx is not a FD object.');
end
if ~strcmp(class(fdobjy), 'fd')
    error('FDOBJx is not a FD object.');
end

%  Center functions if required

if nargin < 6, centerfns = 1; end

if centerfns
    fdobjx = center(fdobjx);
    fdobjy = center(fdobjy);
end

%  check that functions have the same number of replications

coefx  = getcoef(fdobjx);
coefy  = getcoef(fdobjy);
coefdx = size(coefx);
coefdy = size(coefy);
nrepx  = coefdx(2);
nrepy  = coefdy(2);

if (nrepx ~= nrepy)
    error('The numbers of replications are not equal.');
end

%  check that there are more one replication

if nrepx < 2
    error('There is only one replication.');
end

nrep = nrepx;

%  get basis information for data functions

fdbasisx = getbasis(fdobjx);
fdbasisy = getbasis(fdobjy);

%  set up default fdPar objects

if nargin < 5,
    ccafdPary = fdPar(fdbasisy, 2, 1e-10);
end
if nargin < 4,
    ccafdParx = fdPar(fdbasisx, 2, 1e-10);
end

%  check fdPar objects

if strcmp(class(ccafdParx), 'fd') || strcmp(class(ccafdParx), 'basis')
    ccafdParx = fdPar(ccafdParx);
end
if strcmp(class(ccafdPary), 'fd') || strcmp(class(ccafdPary), 'basis')
    ccafdPary = fdPar(ccafdPary);
end

if ~strcmp(class(ccafdParx), 'fdPar')
    error('CCFDPAROBJx is not an FDPAR object.')
end
if ~strcmp(class(ccafdPary), 'fdPar')
    error('CCFDPAROBJy is not an FDPAR object.')
end

%  get basis objects for weight functions

wtbasisx = getbasis(getfd(ccafdParx));
wtbasisy = getbasis(getfd(ccafdPary));

%   Set up essential cross product matrices

Jmatx = inprod_basis(fdbasisx, wtbasisx);
Jmaty = inprod_basis(fdbasisy, wtbasisy);

Jx    = coefx'*Jmatx;
Jy    = coefy'*Jmaty;
PVxx  = Jx' * Jx./nrep; 
PVyy  = Jy' * Jy./nrep;

%  add penalty if either LAMBDA positive

lambdax = getlambda(ccafdParx);
lambday = getlambda(ccafdPary);

if lambdax > 0
    Lfdobjx = getLfd(ccafdParx);
    Kmatx   = eval_penalty(wtbasisx, Lfdobjx);
    PVxx    = PVxx + lambdax * Kmatx;
end
if lambday > 0 
    Lfdobjy = getLfd(ccafdPary);
    Kmaty   = eval_penalty(wtbasisy, Lfdobjy);
    PVyy    = PVyy + lambday * Kmaty;
end

%  set up matrix to be analyzed

Vxy   = Jx' * Jy./nrep;

%  do eigenanalysis

geigstr = geigen(Vxy, PVxx, PVyy);

%  set up canonical correlations and coefficients for weight functions

if nargin < 3, ncan = 2; end

canwtcoefx = geigstr.Lmat(:,1:ncan);
canwtcoefy = geigstr.Mmat(:,1:ncan);

corrs = diag(geigstr.values);

%   Normalize the weight functions

for j = 1:ncan
    temp = squeeze(canwtcoefx(:,j));
    temp = temp./sqrt(sum(temp.^2));
    canwtcoefx(:,j) = temp;
    temp = squeeze(canwtcoefy(:,j));
    temp = temp./sqrt(sum(temp.^2));
    canwtcoefy(:,j) = temp;
end

%  set up final results in struct object CCASTR

wtnamesx  = getnames(fdobjx);
wtnamesy  = getnames(fdobjy);
wtnamesx{3} = ['Can. wt. fn. for ',wtnamesx{3}];
wtnamesy{3} = ['Can. wt. fn. for ',wtnamesy{3}];
ccawtfdx  = fd(canwtcoefx, wtbasisx, wtnamesx);
ccawtfdy  = fd(canwtcoefy, wtbasisy, wtnamesy);

ccavarx = Jx * canwtcoefx;
ccavary = Jy * canwtcoefy;

ccastr.wtfdx = ccawtfdx;
ccastr.wtfdy = ccawtfdy;
ccastr.varx  = ccavarx;
ccastr.vary  = ccavary;
ccastr.corrs = corrs;
