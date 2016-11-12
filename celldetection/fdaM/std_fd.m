function stdfd = std_fd(fdobj)
%  STD  Compute standard deviation functions for functional observations.
%  Argument:
%  FD    ... a functional data object
%  Return:
%  STDFD ... a functional data for the standard deviation functions

%  This version is identical to that defined in std.m, and is
%  included for compatibility with the R function std.fd

%  last modified 27 July 2012

  coef     = getcoef(fdobj);
  coefd    = size(coef);
  ndim     = length(coefd);
  if (coefd(1) == 1)
    error('Only one replication found.');
  end

  basisfd  = getbasis(fdobj);
  nbasis   = getnbasis(basisfd);
  rangeval = getbasisrange(basisfd);

  varbifd  = var(fdobj);

  neval    = 10*nbasis + 1;
  evalarg  = linspace(rangeval(1), rangeval(2), neval)';
  vararray = eval_bifd(varbifd, evalarg, evalarg);
  nvdim    = length(size(vararray));

  if (ndim == 2)
    stdmat  = sqrt(diag(vararray));
  else
    nvar = coefd(3);
    stdmat = zeros(neval, nvar);
    m = 0;
    for j = 1:nvar
      m = m + j;
      if (nvdim == 3)
        stdmat(:,j) = sqrt(diag(vararray(:,:,1,m)));
      else
        stdmat(:,j) = sqrt(diag(vararray(:,:,m)));
      end
    end
  end
  stdcoef = project_basis(stdmat, evalarg, basisfd);

  fdnames = getnames(fdobj);
  fdnames{2} = 'Std. Dev.';
  fdnames{3} = ['Std. Dev. ',fdnames{3}];

  stdfd = fd(stdcoef, basisfd, fdnames);


