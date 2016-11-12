function isadouble = isa_double(x)
%  ISA_DOUBLE checks X for having the DOUBLE class

%  last modified 20 September 2009

isadouble = strcmp(class(x), 'double');