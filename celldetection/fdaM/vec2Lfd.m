function Lfdobj = vec2Lfd(bwtvec, rangeval)
%VEC2LFD converts a vector of length m to a linear differential
%  operator object of order m.  The range of the
%  functional data object in any cell is set to RANGEVAL.  
%  In the event that BWTVEC is already a linear differential operator
%  object, it returns the object.  

%  Last modified 30 September 2009

%  return BWTVEC if it is of class LFD

if isa_Lfd(bwtvec) 
    Lfdobj = bwtvec;
    return
end

%  check BWTVEC

if ~isnumeric(bwtvec) 
    error('Argument not a vector and not a linear differential operator.');
end
 
m = length(bwtvec);

%  set default range

if nargin < 2,  rangeval = [0,1];  end

%  set up the list object for the homogeneous part

if m == 0
    %  if derivative is zero, BWTLIST is empty
    bwtcell = [];
else 
    conbasis = create_constant_basis(rangeval);
    bwtcell  = cell(m,1);
    for j = 1:m 
        bwtcell{j} = fdPar(fd(bwtvec(j), conbasis));
    end
end

%  define the Lfd object

Lfdobj = Lfd(m, bwtcell);


