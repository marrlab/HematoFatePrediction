function isafd = isa_fd(x)
%  ISA_FD checks X for having the FD class

%  last modified 20 September 2009

isafd = strcmp(class(x), 'fd');

