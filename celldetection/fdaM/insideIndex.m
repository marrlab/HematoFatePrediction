function ind = insideIndex(X, Y, p, t, tricoef)
%  insideIndex returns the index of the triangle containing the point
% (X,Y) if such a triangle exists, and NaN otherwise.
%  TRICOEF may have already been calculated for efficiency,
%  but if the function is called with four arguments, it is calculated.

small = 10000*eps;

ntri   = size(t,1);
indtri = (1:ntri)';

%  compute coefficients for computing barycentric coordinates if
%  needed

if nargin < 5
    tricoef = zeros(ntri,4);
    tricoef(:,1) = p(t(:,1),1)-p(t(:,3),1);
    tricoef(:,2) = p(t(:,2),1)-p(t(:,3),1);
    tricoef(:,3) = p(t(:,1),2)-p(t(:,3),2);
    tricoef(:,4) = p(t(:,2),2)-p(t(:,3),2);
    detT = tricoef(:,1).*tricoef(:,4) - tricoef(:,2).*tricoef(:,3);
    tricoef = tricoef./(detT*ones(1,4));
end

%  compute barycentric coordinates

r3 = X - p(t(:,3),1);
s3 = Y - p(t(:,3),2);
lam1 = ( tricoef(:,4).*r3 - tricoef(:,2).*s3);
lam2 = (-tricoef(:,3).*r3 + tricoef(:,1).*s3);
lam3 = 1 - lam1 - lam2;

%  test these coordinates for a triple that are all between 0 and 1

int  = (-small <= lam1 & lam1 <= 1+small) & ...
       (-small <= lam2 & lam2 <= 1+small) & ...
       (-small <= lam3 & lam3 <= 1+small);

%  return the index of this triple, or NaN if it doesn't exist

indi = indtri(int);
if isempty(indi)
    ind = NaN;
else
    ind = min(indi);
end
