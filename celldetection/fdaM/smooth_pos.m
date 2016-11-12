function [Wfd, Fstr] = ...
    smooth_pos(argvals, y, fdParobj, wtvec, conv, iterlim, dbglev)
%SMOOTH_POS smooths the relationship of Y to ARGVALS by fitting 
%     a positive function of the form
%                   f(t) = exp W(t)
%     where  W  is a function defined over the same range as ARGVALS,
%                     W  = log f.
%  The fitting criterion is penalized mean squared error:
%    PENSSE(lambda) = \sum w_i[y_i - f(t_i)]^2 +
%                     \lambda * \int [L W]^2 
%  The function W(t) is expanded by the basis in functional data object
%    Wfd.   The coefficients of this expansion are called "coefficients"
%    in the comments, while the b's are called "regression coefficients"
%
%  Arguments are:
%  ARGVALS ...  Argument value array of length N, where N is the number of 
%               observed curve values for each curve.  It is assumed that
%               that these argument values are common to all observed
%               curves.  If this is not the case, you will need to 
%               run this function inside one or more loops, smoothing
%               each curve separately.
%  Y       ...  Function value array (the values to be fit).
%               If the functional data are univariate, this array will be 
%               an N by NCURVE matrix, where N is the number of observed
%               curve values for each curve and NCURVE is the number of
%               curves observed.
%               If the functional data are muliivariate, this array will be 
%               an N by NCURVE by NVAR matrix, where NVAR the number of
%               functions observed per case.  For example, for the gait
%               data, NVAR = 2, since we observe knee and hip angles.
%  FDPAROBJ ... A functional parameter or fdPar object.  This object 
%               contains the specifications for the functional data
%               object to be estimated by smoothing the data.  See
%               comment lines in function fdPar for details.
%               The functional data object WFD in FDPAROBJ is used
%               to initialize the optimization process.
%               It's coefficient array contains the starting values for 
%               the iterative minimization of mean squared error.
%  WTVEC   ...  a vector of weights, a vector of N one's by default.
%  CONV    ... convergence criterion
%  ITERLIM ...  maximum number of iterations, 50 by default.
%  DBGLEV  ...  Controls the level of output on each iteration.  If 0,
%               no output, if 1, output at each iteration, if higher, 
%               output at each line search iteration. 1 by default.
%
%  Returns:
%  WFD     ...  Functional data object for W. 
%               Its coefficient matrix an N by NCURVE (by NVAR) matrix
%               (or array), depending on whether the functional
%               observations are univariate or multivariate.
%  FSTR ...  A struct object or a cell array of struct objects, one for 
%            each curve (and each variable if functions are multivariate).
%            Each struct object has slots:
%                 f    ... The sum of squared errors   
%                 grad ... The gradient  
%                 norm ... The norm of the gradient  

%  last modified 21 September 2008

if nargin < 3
    error('The first three arguments are not specified.');
end

%  check ARGVALS

[argvals, n] = argcheck(argvals);

%  check Y and compute root-mean-square of values

[y, ncurve, nvar, ndim] = ycheck(y, n);
 
%  check fdParobj

fdParobj = fdParcheck(fdParobj);
lambda   = getlambda(fdParobj);

%  set up LFDOBJ

Lfdobj = getLfd(fdParobj);
Lfdobj = int2Lfd(Lfdobj);

%  set up BASIS

Wfd0     = getfd(fdParobj);
basisobj = getbasis(Wfd0);
nbasis   = getnbasis(basisobj);
active   = 1:nbasis;

%  set up some intial arrays

coef0 = getcoef(Wfd0);  %  initial coefficients

%  set some default arguments and constants

if nargin < 7, dbglev  = 1;          end
if nargin < 6, iterlim = 50;         end
if nargin < 5, conv    = 1e-4;       end
if nargin < 4, wtvec   = ones(n,1);  end
if nargin < 3
    error('Less than three arguments are supplied.');
end

%  check WTVEC

wtvec = wtcheck(n, wtvec);

%  set up some arrays

climit = [-50,0;0,400]*ones(2,nbasis);
dbgwrd = dbglev > 1;

%  initialize matrix Kmat defining penalty term

if lambda > 0
    Kmat = lambda.*eval_penalty(basisobj, Lfdobj);
else
    Kmat = zeros(nbasis);
end

%  --------------------------------------------------------------------
%              loop through variables and curves
%  --------------------------------------------------------------------

%  set up arrays and cell arrays to contain returned information

if ndim == 2
    coef = zeros(nbasis,ncurve);
else
    coef = zeros(nbasis,ncurve,nvar);
end

if nargout > 2
    if ncurve > 1 || nvar > 1
        Fstr = cell(ncurve,nvar);
    else
        Fstr = [];
    end
end

for ivar=1:nvar
    if ndim == 2
        sclfac = mean(y(:).^2);
    else
        sclfac = mean(mean(y(:,:,ivar)).^2);
    end
    for icurve=1:ncurve
        if ndim == 2
            yi    = squeeze(y(:,icurve));
            cveci = squeeze(coef0(:,icurve));
        else
            yi    = squeeze(y(:,icurve,ivar));
            cveci = squeeze(coef0(:,icurve,ivar));
        end

        %  evaluate log likelihood
        %    and its derivatives with respect to these coefficients

        [PENSSE, DPENSSE] = PENSSEfun(argvals, yi, wtvec, basisobj, ...
                                      cveci, Kmat);

        %  compute initial badness of fit measures

        gvec         = DPENSSE;
        Foldstr.f    = PENSSE;
        Foldstr.grad = gvec;
        Foldstr.norm = sqrt(mean(gvec.^2));

        %  compute the initial expected Hessian

        hmat = PENSSEhess(argvals, yi, wtvec, basisobj, cveci, Kmat);

        %  evaluate the initial update vector for correcting initial bmat

        deltac = -hmat\gvec;

        %  initialize iteration status arrays

        iternum = 0;
        status = [iternum, Foldstr.f, Foldstr.norm];
        if dbglev >= 1
            fprintf(['\nResults for curve ',num2str(icurve), ...
                ' and variable ',num2str(ivar),'\n'])
            fprintf('\nIteration  Criterion  Grad. Norm\n')
            fprintf('\n%5.f     %10.4f %10.4f\n', status);
        end
        if dbglev > 2
            for ibasis = 1:nbasis, fprintf('%10.4f%', cveci(ibasis)); end
            fprintf('\n');
        end

        %  ---------------------  Begin main iterations  ---------------

        STEPMAX = 5;
        MAXSTEP = 400;
        trial   = 1;
        linemat = zeros(3,5);

        %  ---------------  beginning of optimization loop  -----------

        for iter = 1:iterlim
            iternum = iternum + 1;
            %  take optimal stepsize
            dblwrd = [0,0]; limwrd = [0,0]; ind = 0;
            %  compute slope
            Fstri = Foldstr;
            linemat(2,1) = sum(deltac.*gvec);
            %  normalize search direction vector
            sdg     = sqrt(sum(deltac.^2));
            deltac  = deltac./sdg;
            linemat(2,1) = linemat(2,1)/sdg;
            %  return with error condition if initial slope is nonnegative
            if linemat(2,1) >= 0
                disp('Initial slope nonnegative.')
                break;
            end
            %  return successfully if initial slope is very small
            if linemat(2,1) >= -sclfac*1e-5;
                if dbglev>1, disp('Initial slope too small'); end
                break;
            end
            linemat(1,1:4) = 0;
            linemat(2,1:4) = linemat(2,1);
            linemat(3,1:4) = Foldstr.f;
            stepiter  = 0;
            if dbglev > 1
                fprintf('      %3.f %10.4f %10.4f %10.4f\n', ...
                    [stepiter, linemat(:,1)']);
            end
            ips = 0;
            %  first step set to trial
            linemat(1,5) = trial;
            %  Main iteration loop for linesrch
            for stepiter = 1:STEPMAX
                %  ensure that step does not go beyond limits on parameters
                %  check the step size
                [linemat(1,5),ind,limwrd] = ...
                    stepchk(linemat(1,5), cveci, deltac, ...
                            limwrd, ind, climit, active, dbgwrd);
                if linemat(1,5) <= 1e-5
                    %  Current step size too small ... terminate
                    Fstri   = Foldstr;
                    cvecnew = cveci;
                    gvecnew = gvec;
                    if dbglev > 1
                        fprintf(  ...
                            'Stepsize too small:  %10.4f\n', linemat(1,5));
                    end
                    break;
                end
                cvecnew = cveci + linemat(1,5).*deltac;
                %  compute new function value and gradient
                [PENSSE, DPENSSE] = ...
                    PENSSEfun(argvals, yi, wtvec, basisobj, cvecnew, Kmat);
                gvecnew    = DPENSSE;
                Fstri.f    = PENSSE;
                Fstri.grad = gvecnew;
                Fstri.norm = sqrt(mean(gvecnew.^2));
                linemat(3,5) = Fstri.f;
                %  compute new directional derivative
                linemat(2,5) = sum(deltac.*gvecnew);
                if dbglev > 1
                    fprintf('      %3.f %10.4f %10.4f %10.4f\n', ...
                        [stepiter, linemat(:,5)']);
                end
                %  compute next step
                [linemat,ips,ind,dblwrd] = ...
                    stepit(linemat, ips, dblwrd, MAXSTEP);
                trial  = linemat(1,5);
                %  ind == 0 implies convergence
                if ind == 0 || ind == 5, break; end
                %  end iteration loop
            end

            cveci = cvecnew;
            gvec  = gvecnew;
            status = [iternum, Fstri.f, Fstri.norm];
            fprintf('%5.f     %10.4f %10.4f\n', status);
            %  test for convergence
            if abs(Fstri.f-Foldstr.f) < conv*sclfac;
                break;
            end
            if Fstri.f >= Foldstr.f, break; end
            %  compute the Hessian
            hmat = PENSSEhess(argvals, yi, wtvec, basisobj, cveci, Kmat);
            %  evaluate the update vector
            deltac   = -hmat\gvec;
            cosangle = -gvec'*deltac/sqrt(sum(gvec.^2)*sum(deltac.^2));
            if cosangle < 0
                if dbglev > 1, disp('cos(angle) negative'); end
                deltac = -gvec;
            end
            Foldstr = Fstri;
        end

        %  save coefficients in arrays COEF
        if ndim == 2
            coef(:,icurve) = cveci;
        else
            coef(:,icurve,ivar) = cveci;
        end
        
         %  save Fstr if required in cell array.
        
        if nargout > 2
            if ncurve == 1 && nvar == 1
                Fstr = Fstri;
            else
                Fstr{icurve,ivar} = Fstri;
            end
        end
        
    end
end

Wfd = fd(coef, basisobj);

%  ---------------------------------------------------------------

function [PENSSE, DPENSSE] = ...
    PENSSEfun(argvals, yi, wtvec, basisobj, cveci, Kmat)
n       = length(argvals);
phimat  = getbasismatrix(argvals, basisobj);
Wvec    = phimat*cveci;
EWvec   = exp(Wvec);
res     = yi - EWvec;
PENSSE  = mean(wtvec.*res.^2) + cveci'*Kmat*cveci;
DPENSSE = -2.*phimat'*(wtvec.*res.*EWvec)./n + 2.*Kmat*cveci;

%  ---------------------------------------------------------------

function D2PENSSE = PENSSEhess(argvals, yi, wtvec, basisobj, cveci, Kmat)
n = length(argvals);
nbasis   = getnbasis(basisobj);
phimat   = getbasismatrix(argvals, basisobj);
Wvec     = phimat*cveci;
EWvec    = exp(Wvec);
res      = yi - EWvec;
Dres     = ((res.*EWvec)*ones(1,nbasis)) .* phimat;
D2PENSSE = 2.*Dres'*((wtvec*ones(1,nbasis)).*Dres)./n + 2.*Kmat;

