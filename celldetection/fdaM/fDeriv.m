function Dfdobj = fDeriv(fdobj, nderiv)
coefmat    = getcoef(fdobj);
[nbasis,N] = size(coefmat);
basisobj   = getbasis(fdobj);
type       = getbasistype(basisobj);
if ~strcmp(type, 'bspline') 
    error('Basis is not of spline type.');
end
rng        = getbasisrange(basisobj);
pars       = getbasispar(basisobj);
nbasis     = getnbasis(basisobj);
norder     = nbasis - length(pars);

knots = [rng(1)*ones(1,norder),pars,rng(2)*ones(1,norder)]';

for ideriv=1:nderiv
    coefmat = diff(coefmat);
    Dknot   = knots(norder+1:nbasis+norder-1) - knots(2:nbasis);
    nbasis  = nbasis-1;
    norder  = norder-1;
    knots   = knots(2:length(knots)-1);
    coefmat = norder.*coefmat./repmat(Dknot,1,N);
end

Dbreaks   = [rng(1),pars,rng(2)];
Dbasisobj = create_bspline_basis(rng,nbasis,norder,Dbreaks);
Dfdobj    = fd(coefmat,Dbasisobj);
