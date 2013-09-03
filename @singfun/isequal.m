function out = isequal(f, g)
%ISEQUAL   Test if SINGFUN objects F and G are equal.
%   ISEQUAL(F, G) returns TRUE if the SINGFUN objects F and G have the same
%   underlyig SMOOTHPART and the same EXPONENTS. The type of singularity
%   at the ends may be differnet but as long as the EXPONENTS agree within
%   SINGFUN tolerance, they are considered equal.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

tol = singfun.pref.singfun.exponentTol;

out = all(abs(f.exponents - g.exponents) < tol) && ...
    isequal(f.smoothPart, g.smoothPart);

end
