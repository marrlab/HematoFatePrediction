function isafdPar = isa_fdPar(fdParstr)
%  ISA_FDPAR  checks an object to see if it is of the FD class 

%  last modified 1 May 2009

  isafdPar = 1;
  if ~(strcmp(class(fdParstr), 'fdPar'))
    isafdPar = 0;
    return;
  end

