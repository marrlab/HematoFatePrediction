function basisobj = create_IRT3PL_basis(rangeval, nbasis, shift, dropind)
%  CREATE_IRT3PL_BASIS  Creates a monomial basis:, x^i_1, x^i_2, ...
%  combined with an additional basis function log(exp(x)+shift)
%  for use with simple item response functions to model the log-odds
%  function.  The number of monomial terms is NBASIS - 1, and the highest
%  exponent is NBASIS - 2
%  Arguments:
%  RANGEVAL ... an array of length 2 containing the lower and upper
%               boundaries for the rangeval of argument values.  If a
%               single value is input, it must be positive and the lower
%               limit of the range is set to 0.
%  NBASIS    ... number of basis functions.  Defaults to 3.
%  SHIFT     ... Shift parameter in log(exp(x)+shift)
%  DROPIND   ... a set of indices in 1:NBASIS of basis functions to drop
%                when basis objects are arguments.  Default is [];
%  Return:
%  BASIS  ... a functional data basis object of type 'IRT3PL'

%  last modified 19 June 2012

%  default RANGEVAL

if nargin < 1, rangeval = [-3,3];  end

%  check RANGEVAL

if length(rangeval) == 1
    if rangeval <= 0
        error('RANGEVAL a single value that is not positive.');
    end
    rangeval = [-rangeval,rangeval];
end

if rangechk(rangeval) ~= 1
    error('RANGEVAL is not a legitimate range.');
end

%  set up default arguments

if nargin < 2, nbasis = 3;    end
if nargin < 3, shift = 1;     end
if nargin < 4, dropind = [];  end

exponents = 0:(nbasis-2);

%  check DROPIND

if length(dropind) > 0
    if length(dropind) >= nbasis
        error('Too many index values in DROPIND.');
    end
    dropind = sort(dropind);
    if length(dropind) > 1
        if any(diff(dropind)) == 0
            error('Multiple index values in DROPIND.');
        end
    end
    for i=1:length(dropind);
        if dropind(i) < 1 || dropind(i) > nbasis
            error('An index value is out of range.');
        end
    end
end


%  set up basis object

type        = 'IRT3PL';
params      = shift;
quadvals    = [];
values      = {};
basisvalues = {};

basisobj = basis(type, rangeval, nbasis, params, ...
                 dropind, quadvals, values, basisvalues);

