function [mon_grad, mon_hess, mon_cross, basiscell] = ...
        mon_grad_fn(y, x, cvec, lambda, basisobj, Kmat, basiscell)
n      = length(x);
nbasis = size(cvec,1);
Wfd    = fd(cvec, basisobj);
[f,  basiscell]  =   monfn(x, Wfd, basiscell);
[Df, basiscell]  = mongrad(x, Wfd, basiscell);
if nargout > 1
    [D2f, basiscell] = monhess(x, Wfd, basiscell);
end    
xmat   = [ones(n,1), f];
%  compute regression coefs.
b  = xmat\y;
b2 = b(2);
%  update fitted values and residuals
mon_grad = zeros(nbasis,1);
for j=1:nbasis
    Dxmatj = Df(:,j);
    yDx    = y'*Dxmatj*b2;
    xDx    = xmat'*Dxmatj;
    mon_grad(j) = b'*(xDx+xDx')*b - 2*yDx;
end
mon_grad = mon_grad/n;
if lambda > 0
    mon_grad = mon_grad + 2.*Kmat*cvec;
    if nargout > 2
        mon_hess = mon_hess + 2.*Kmat;
    end
end

