function cormat = cor_fd(evalarg1, fdobj1, evalarg2, fdobj2)

%  compute correlation matrix / matrices for functional observations

%  Last modified 6 January 2008 by Jim Ramsay
%
% 1.  Compute var1 = bivariate data object for variance(fdobj1)
% %
var1 = var_fd(fdobj1);
%
% 2.  Evaluate var1 at evalarg1
%
evalVar1 = eval_bifd(evalarg1, evalarg1, var1);
%if(   if
% 3.  If missing(fdobj2) convert evalVar1 to correlations
%
if nargin < 4
    sizeV1 = size(evalVar1);
    ndV1  = length(sizeV1);
    if ndV1 < 3
        s1 = sqrt(diag(evalVar1));
        cormat = evalVar1./(s1*s1');
    else
        if sizeV1(3) ~= 1
            error(['Bug in cor_fd:  Programmed only for ', ...
                'matrices or 4-d arrays with size(3)==1.  oops.'])
        end
        size1 = size(getcoef(fdobj1));
        nVars = size1(3);
        %         The following identifies the levels of evalVar1
        %         containing the covariance matrix of each variable with itself
        evalV_diag = 2:2:(nVars+1);
        %         Compute the standard deviation vector for each variable
        nPts = length(evalarg1);
        s1 = zeros(nPts, nVars);
        for i = 1:nVars
            s1(:,i) = sqrt(diag(squeeze(evalVar1(:,:,1,evalV_diag(i)))));
        end
        %         Now compute the correlations
        cormat = evalVar1;
        m = 0;
        for i = 1:nVars
            for j = 1:i
                m = m+1;
                cormat(:,:,1,m) = evalVar1(:,:,1,m)./(s1(:,i)*s1(:,j))';
            end
            return
        end
    end
else
    %
    % 4.  fdobj2 was also provided
    %
    %  4.1  var_df(fdobj2)
    var2 = var_fd(fdobj2);
    %  4.2.  Evalate at evalarg2
    evalVar2 = eval_bifd(evalarg2, evalarg2, var2);
    %  4.3.  var12 cross covariance
    %*** If fdobj1 or fdobj2 are multivariate, var_fd will compla:.
    var12 = var_fd(fdobj1, fdobj2);
    %  4.4.  Evaluate the cross covariances
    evalVar12 = eval_bifd(evalarg1, evalarg2, var12);
    %  4.5.  Convert evalVar12 to correlations
    s1 = sqrt(diag(evalVar1));
    s2 = sqrt(diag(evalVar2));
    cormat = evalVar12./(s1*s2');
end

