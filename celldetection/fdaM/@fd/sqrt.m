function fdaroot = sqrt(fdobj)
%  Compute the square root of a functional data object.

%  Last modified 27 October 2009

if ~isa_fd(fdobj)
    error('Argument is not a functional data object.');
end

fdaroot = fdobj.^0.5;