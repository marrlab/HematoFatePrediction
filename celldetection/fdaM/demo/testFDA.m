%% Introducing the FDA Toolbox
%
% Disclaimer -- the FDA toolbox is not an official Matlab toolbox; it is a
% collection of privately-written functions and comes with no garuantees. 
%
% To get the toolbox
%
% 1. Visit functionaldata.org
%
% 2. Click on software and follow the link at the bottom of the page
%
% 3. Alternatively, go straight to
%
%  ftp://ego.psych.mcgill.ca/pub/ramsay/FDAfuns/Matlab
%
% 4. Dowload     Matlabfunctions.zip
%
% 5. Follow the instruction in INSTALL.MATLAB.WIN.txt
%
%     Linux users; the installation is analagous
%
%     I keep my fda package in
%
%       F://work/External/matlab_packages/fdaM
%
% Now in order to make use of these functions, you need to add them to your
% path:


% Now let's load some data

load 'handwriting.mat'

% Now we would like to start smoothing the data

%% Defining a Basis

help basis

% All basis objects require a range to be specified:

basis_range = [min(fdatime) max(fdatime)];

% Starting with a B-spline basis:

breaks = linspace(basis_range(1),basis_range(2),11);

figure(1)
plot(breaks,'.')

% I'll choose cubic B-splines

norder = 4;

% And that gives 

nbasis = length(breaks) + norder - 2;

% This is enough to define the basis

bspline_obj = create_bspline_basis(basis_range,nbasis,norder,breaks);

% Now that I've got a basis, I can evaluate it:

basisvals = eval_basis(fdatime,bspline_obj);

size(basisvals)

figure(2)
plot(fdatime,basisvals)

% There is also a shortcut

figure(3)
plot(bspline_obj)


% Lets look at just one

figure(4)
plot(fdatime,basisvals(:,4),'linewidth',2)
hold on
plot([breaks; breaks],[-ones(1,length(breaks)); ones(1,length(breaks))],'r--') 


% Now I can evaluate some derivatives; the third argument in eval_basis
% gives an order of derivative to evaluate. 

bbasisvals = eval_basis(fdatime,bspline_obj,1);
b2basisvals = eval_basis(fdatime,bspline_obj,2);

plot(fdatime,bbasisvals(:,7)/3,'m','linewidth',2)
plot(fdatime,b2basisvals(:,10)/40,'k','linewidth',2)
hold off

% Here I have had to scale the derivatives in order to make them comparable
% in size. This is important numerically -- using a very fine basis can
% mean that taking derivatives leads to numerical overflow or underflow. 


% To give a comparison, let's try a linear basis

norder = 2;
nbasis = length(breaks);

bspline_obj = create_bspline_basis(basis_range,nbasis,norder,breaks);
basisvals = eval_basis(fdatime,bspline_obj);

figure(5)
plot(fdatime,basisvals)

axis([0 2.3 0 2])

% There is probably not enough resolving power in this basis, so I'll try
% more knots

norder = 4;

breaks = linspace(basis_range(1),basis_range(2),101);
nbasis = length(breaks) + norder - 2;

bspline_obj = create_bspline_basis(basis_range,nbasis,norder,breaks);

basisvals = eval_basis(fdatime,bspline_obj);

figure(6)
plot(fdatime,basisvals)
axis([0 2.3 0 1])

% Lets try a fourier basis:

nbasis_f = 10;
period = basis_range(2)-basis_range(1);

fourier_obj = create_fourier_basis(basis_range,nbasis_f,period);

fbasisvals = eval_basis(fdatime,fourier_obj);

figure(7)
plot(fdatime,fbasisvals(:,1:5));

% And maybe a polynomial basis on a centered range

nbasis_p = 10;
exponents = 0:9;    % Note that I do not need to use integer exponents here

monomial_obj = create_monomial_basis(basis_range,nbasis_p,exponents);

mbasisvals = eval_basis(fdatime,monomial_obj);

figure(8)
plot(fdatime,mbasisvals)


% Accessing aspects of a basis object

getbasisrange(bspline_obj)   % range of the basis

getnbasis(bspline_obj)       % number of basis functions

bvals = getbasismatrix(fdatime,fourier_obj);   % same as eval_basis
size(bvals)

%% Creating Functional Data Objects
%
% fd objects combine a basis with coefficients

coefs = randn(getnbasis(fourier_obj),1);

fd_obj = fd(coefs,fourier_obj);

figure(9)
plot(fd_obj)

% Alternatively

fvals = eval_fd(fdatime,fd_obj);

figure(10)
plot(fdatime,fvals);


% We can also define collections of fd objects

coefs = randn(getnbasis(fourier_obj),10);

fdnames = {'fdatime','reps','height'};

fd_obj = fd(coefs,fourier_obj,fdnames);

figure(11)
plot(fd_obj)

% I can evaluate derivatives by

fvals = eval_fd(fdatime,fd_obj,3);

size(fvals)

figure(12)
plot(fdatime,fvals);


%% Smoothing Data

% We saw in class that smoothing with a basis expansion just requires
% solving the normal equations. There is a function to do this
% automatically:

fdax_b = data2fd(fdarray(:,:,1),fdatime',bspline_obj);

figure(13)
plot(fdax_b)

figure(12)
plotfit_fd(fdarray(:,:,1),fdatime',fdax_b)

% Now if I want to plot a derivative:

figure(14)
plot(fdax_b,1)

% compare this with

figure(15)
plot(fdatime(1:1400),diff(fdarray(:,:,1)))

% Let's look at a fourrier basis:

fdax_f = data2fd(fdarray(:,:,1),fdatime,fourier_obj);

plotfit_fd(fdarray(:,:,1),fdatime,fdax_f);


% Or a monomial basis

fdax_m = data2fd(fdarray(:,:,1),fdatime,monomial_obj);

plotfit_fd(fdarray(:,:,1),fdatime,fdax_m);


% I can recover various

coef_f = getcoef(fdax_f);

basis_f = getbasis(fdax_f);

basis_f

%% Numerical Quadrature

% For approximating quadrature; usually. 

% I can set quadrature points and weights. 

r_b = getbasisrange(basis_f);

m = 1000;

quadpts = linspace(r_b(1),r_b(2),2*m+1);  % Define a fine grid over the range

quadwts = ones(2*m+1,1);                  % Simpson's rule quadrature weights
quadwts(2*(1:(m-1))) = 2;
quadwts(2*(1:m)+1) = 4;

quadwts = quadwts*(r_b(2)-r_b(1))/(6*m);

basis_f = putquadvals(basis_f,[quadpts' quadwts]);

quads = getquadvals(basis_f);

size(quads)

% There are times when it may be useful to have quick access to derivatives
% of the basis values at the quadrature points.

k = 2;  % number of derivatives I'm interested in

values = cell(3,1);

for i = 0:k, values{i+1} = diag(sqrt(quadwts))*eval_basis(quadpts,basis_f,i); end

basis_f = putvalues(basis_f,values);

vals = getvalues(basis_f,1);

size(vals)

figure(16)
plot(vals)

% Note that I multiplied by the square-root of quadwts. This is because I
% will usually be interested in estimating the integral of one basis
% against another.
%
% If I have B1 and B2 evaluated at quadrature points, then their inner
% product is 
%
% B1'*diag(quadwts)*B2
%
% multiplying by the square root of quadwts avoids using the middle term. 

%% Multiple Knots

% In order to demonstrate the use of B-splines for known discontinuities

fda = [fdarrayx; fdarrayy];

fdatime2 = [fdatime; max(fdatime)+ fdatime(2) + fdatime];

% Now set up a order-4 B-splines with three knots at 2.3. 

r_f = [min(fdatime2) max(fdatime2)];

breaks = linspace(r_f(1),r_f(2),201);

knots = sort([breaks 2.3*ones(1,2)]);

norder = 4;

nbasis = length(knots)+norder-2;

basis_b = create_bspline_basis(r_f,nbasis,norder,knots);

fd_obj = data2fd(fda(:,1),fdatime2,basis_b);

fda1 = eval_fd(fdatime2,fd_obj);

figure(17)
plot(fdatime2,fda1)



