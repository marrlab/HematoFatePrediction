function [grd, Hess] = grdHessm(f,x,toler) 
%Input
%       - f is the function
%       - x is the differentiation point (vector valued)
%       - toler is the tolerance for the relative error
%Output 
%       - grd is the gradient
%       - Hess is the Hessian

%  last modified 13 September 2009 by Jim

%  set default for toler

if nargin < 3, toler = 1e-4;  end

p = length(x);

grd   = zeros(p,1);
dHess = zeros(p,1);
Hess  = zeros(p,p);
for i = 1:p
    [grdi, dHessi] = grdHessD(f,x,i,toler);
    grd(i)   = grdi;
    dHess(i) = dHessi;
end
idx = [kron((1:p),ones(1,p))',repmat((1:p),1,p)'];
idx = idx((idx(:,1)-idx(:,2)~=0),:);
for i = 1:p
    Hess(idx(i,1),idx(i,2)) = HessOD(f,x,idx(i,1),idx(i,2),toler);
    Hess(idx(i,2),idx(i,1)) = Hess(idx(i,1),idx(i,2));
end
for i = 1:p
    Hess(i,i) = dHess(i);
end

function [grdi, dHessi] = grdHessD(f,x,idx,toler)
machine_eps = 2.220446e-16;
p     = length(x);
rerr1 = 1;
rerr2 = 1;
h     = 1;
j     = 1;
hp    = zeros(p,1);
hp(idx) = h;
Dmat1 = zeros(12);
Dmat2 = Dmat1;
fh = f(x+hp);
ff = f(x);
fmh = f(x-hp);
Dmat1(1,1) = (fh-fmh)/(2*h);
Dmat2(1,1) = (fh+2*ff-fmh)/(h^2);
while ((rerr1 > toler) || (rerr2 > toler)) && (j < 12)
    h = h/2;
    hp(idx) = h;
    fh = f(x+hp);
    ff = f(x);
    fmh = f(x-hp);
    Dmat1(j+1,1) = (fh-fmh)/(2*h);
    Dmat2(j+1,1) = (fh-2*ff+fmh)/(h^2);
    for k = 1:j
        Dmat1(j+1,k+1) = Dmat1(j+1,k)+(Dmat1(j+1,k)-Dmat1(j,k))/((4^k)-1);
        Dmat2(j+1,k+1) = Dmat2(j+1,k)+(Dmat2(j+1,k)-Dmat2(j,k))/((4^k)-1);
    end
    err1  = abs(Dmat1(j+1,j+1) - Dmat1(j,j));
    err2  = abs(Dmat2(j+1,j+1) - Dmat2(j,j));
    rerr1 = 2*err1/(abs(Dmat1(j+1,j+1)) + abs(Dmat1(j,j)) + machine_eps);
    rerr2 = 2*err2/(abs(Dmat2(j+1,j+1)) + abs(Dmat2(j,j)) + machine_eps);
    j = j+1;
end
grdi   = Dmat1(j,j);
dHessi = Dmat2(j,j);

function Hessik = HessOD(f,x,idx1,idx2,toler)
machine_eps = 2.220446e-16;
p    = length(x);
rerr = 1;
h    = 1;
j    = 1;
hp1  = zeros(p,1);
hp2  = hp1;
hp1(idx1) = h;
hp2(idx2) = h;
Dmat = zeros(12);
Dmat(1,1) = ...
    (f(x+hp1+hp2)+f(x-hp1-hp2)-f(x-hp1+hp2)-f(x+hp1-hp2))/(4*h^2);
while (rerr > toler) && (j < 12)
    h = h/2;
    hp1(idx1) = h;
    hp2(idx2) = h;
    Dmat(j+1,1) = ...
        (f(x+hp1+hp2)+f(x-hp1-hp2)-f(x-hp1+hp2)-f(x+hp1-hp2))/(4*h^2);
    for k = 1:j
        Dmat(j+1,k+1) = Dmat(j+1,k)+(Dmat(j+1,k)-Dmat(j,k))/((4^k)-1);
    end
    err = abs(Dmat(j+1,j+1)-Dmat(j,j));
    rerr = 2*err/(abs(Dmat(j+1,j+1))+abs(Dmat(j,j))+machine_eps);
    j = j+1;
end
Hessik = Dmat(j,j);


