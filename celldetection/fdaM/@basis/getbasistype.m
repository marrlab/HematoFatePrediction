function typestr = getbasistype (basisobj)
%  GETBASISTYPE   Extracts the type of basis from basis object BASISOBJ.
%    Variations in spelling are possible.

%  last modified 24 May 2013

if ~isa_basis(basisobj)
    error('Argument is not a functional basis object.');
end

switch basisobj.type
    case 'Fourier'
        typestr = 'fourier';
    case 'fourier'
        typestr = 'fourier';
    case 'Fou'
        typestr = 'fourier';
    case 'fou'
        typestr = 'fourier';
    case 'sine'
        typestr = 'sine';
    case 'sin'
        typestr = 'sine';
    case 'Sine'
        typestr = 'sine';
    case 'Sin'
        typestr = 'sine';
    case 'cosine'
        typestr = 'cosine';
    case 'cos'
        typestr = 'cosine';
    case 'Cosine'
        typestr = 'cosine';
    case 'Cos'
        typestr = 'cosine';
    case 'bspline'
        typestr = 'bspline';
    case 'Bspline'
        typestr = 'bspline';
    case 'Bsp'
        typestr = 'bspline';
    case 'bsp'
        typestr = 'bspline';
    case 'monom'
        typestr = 'monom';
    case 'mon'
        typestr = 'monom';
    case 'monomial'
        typestr = 'monom';
    case 'nspline'
        typestr = 'nspline';
    case 'naturalbspline'
        typestr = 'nspline';
    case 'NaturalBspline'
        typestr = 'nspline';
    case 'Nbsp'
        typestr = 'nspline';
    case 'nbsp'
        typestr = 'nspline';
    case 'power'
        typestr = 'power';
    case 'pow'
        typestr = 'power';
    case 'polyg'
        typestr = 'polyg';
    case 'polygon'
        typestr = 'polyg';
    case 'polygonal'
        typestr = 'polyg';
    case 'expon'
        typestr = 'expon';
    case 'exp'
        typestr = 'expon';
    case 'exponential'
        typestr = 'expon';
    case 'const'
        typestr = 'const';
    case 'con'
        typestr = 'const';
    case 'constant'
        typestr = 'const';
    case 'QW'
        typestr = 'QW';
    case 'QWM'
        typestr = 'QWM';
    case 'QS'
        typestr = 'QS';
    case 'slide'
        typestr = 'slide';
    case 'fd'
        typestr = 'fd';
    case 'FEM'
        typestr = 'FEM';
    case 'TP'
        typestr = 'TP';
    case 'fdVariance'
        typestr = 'fdVariance';
    case 'IRT3PL'
        typestr = 'IRT3PL';
    otherwise
        error('Unknown type encountered');
end

