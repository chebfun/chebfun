function out = isequal(f, g)
%ISEQUAL   Test if DELTAFUN objects F and G are equal.
%   ISEQUAL(F, G) returns TRUE if the DELTAFUN objects F and G have the same
%   underlyig FUNPART and the same delta functions.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

tol = singfun.pref.singfun.exponentTol;

out = (f.delta.location - g.delta.location) < 
out = f.funPart == g.funPart;
end
