function regfdobj = synch(time, fdParobj, WfdParobj, hinvwrd, dbglev)
%  Synchronize or register the functions in fdPar object using
%  the W functions in fdPar object WfdParobj over time values in
%  vector time.

%  Arguments:

%  TIME      ...  A vector of N time values
%  FDPAROBJ  ...  Either a fd object of a fdPar object containing
%                 functions to be registered.  It is wise to supply
%                 a complete fdPar object, though, so as to control
%                 the smoothing process in the call to function
%                 smooth_basis below.
%  WFDPAROBJ ...  Either a fd object of a fdPar object containing
%                 the W functions estimated by a registration process.
%                 Again, it is wise to supply
%                 a complete fdPar object, though, so as to control
%                 the smoothing process in the call to function
%                 smooth_morph below.
%  Returns:

%  REGFDOBJ  .... The synchronized or registered functions

%  Last modified 4 August 2008

%  set default level of dbgwrd

if nargin < 5, dbglev = 0;  end

%  set default value of hinvwrd

if nargin < 4, hinvwrd = 0;  end

%  check time

if strcmp(class(time), 'double')
    if ~(size(time,1) == 1 || size(time,2) == 1)
        error('Argument time is not a vector.');
    else
        time = time(:);
    end
else
    error('Argument time is not a vector.');
end

%  check fdParobj

if ~strcmp(class(fdParobj), 'fdPar')
    if strcmp(class(fdParobj), 'fd')
        fdParobj = fdPar(fdParobj);
    else
        error(['Argument fdParobj is neither a fdPar object ', ...
               'nor a fd object.']);
    end
end

%  check WfdParobj

if ~strcmp(class(WfdParobj), 'fdPar')
    if strcmp(class(WfdParobj), 'fd')
        WfdParobj = fdPar(WfdParobj);
    else
        error(['Argument WfdParobj is neither a fdPar object ', ...
               'nor a fd object.']);
    end
end

% set up needed variables

fdobj   = getfd(fdParobj);
coefmat = getcoef(fdobj);
fdmat   = eval_fd(time, fdobj);
if length(size(coefmat)) == 3
    nvar = size(coefmat,3);
else
    nvar = 1;
end
N = size(coefmat,2);
n = length(time);
range = getbasisrange(getbasis(fdobj));
T0 = range(1);
T1 = range(2);

%  extract functions W(t) defining time warping

Wfdobj  = getfd(WfdParobj);

%  loop through curves to be synchronized

coefregmat = coefmat;
for i=1:N
    if dbglev >= 1 && N > 1
        fprintf(['\n\nCurve ',num2str(i),'\n'])
    end;
    %  compute the warping function h(t)
    if hinvwrd
        hi = eval_fd(time, Wfdobj(i));
    else
        hi = monfn(time, Wfdobj(i));
    end
    hi = T0 + (T1-T0).*hi./hi(n);
    hi(1)   = T0;
    hi(n)   = T1;
    %  compute the functional inverse of h(t)
    WfdPari = putfd(WfdParobj, Wfdobj(i));
    Wfdinv  = smooth_morph(hi, time, WfdPari);
    hinv    = monfn(time, Wfdinv);
    hinv    = T0 + (T1-T0).*hinv./hinv(n);
    hinv(n) = T1;
    hinv(1) = T0;
    %  smooth the data over vector hinv to get the registered
    %  functions
    if nvar == 1
        fdmati = squeeze(fdmat(:,i));
        regfdi = smooth_basis(hinv, fdmati, fdParobj);
        coefregmat(:,i) = getcoef(regfdi);
    else
        fdmati = squeeze(fdmat(:,i,:));
        regfdi = smooth_basis(hinv, fdmati, fdParobj);
        coefregmat(:,i,:) = squeeze(getcoef(regfdi));
    end
end

%  set up the registered functions

regfdobj = putcoef(fdobj, coefregmat);