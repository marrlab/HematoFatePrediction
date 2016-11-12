function isabifdPar = isa_bifdPar(bifdParstr)
%  ISA_FDPAR  checks an object to see if it is of the FD class 

%  last modified 1 May 200

  isabifdPar = 1;
  if ~(strcmp(class(bifdParstr), 'bifdPar'))
    isabifdPar = 0;
    return;
  end

