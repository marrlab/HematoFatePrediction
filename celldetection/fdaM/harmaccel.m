function harmaccelLfd = harmaccel(basisobj)
%  Creates a a linear differential operator object of class Lfd.
%  The basis must be a Fourier basis.

if ~strcmp(class(basisobj), 'basis')
    error('BASISOBJ is not a basis object.');
end

type = getbasistype(basisobj);

if ~strcmp(type, 'fourier')
    error('BASISOBJ is not a Fourier basis object.');
end

period   = getbasispar(basisobj);
rangeval = getbasisrange(basisobj);

Lbasisobj    = create_constant_basis(rangeval);
Lcoef        = [0, (2*pi/period)^2, 0];
wfd          = fd(Lcoef, Lbasisobj);
wfdcell      = fd2cell(wfd);     % convert the FD object to a cell object
harmaccelLfd = Lfd(3, wfdcell);  %  define the operator object
