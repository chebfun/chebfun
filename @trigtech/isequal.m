function out = isequal(f, g)
%ISEQUAL   Test if TRIGTECH objects are equal.
%   ISEQUAL(F, G) returns TRUE if the TRIGTECH objects F and G have the same
%   length, values, and coefficients. They may have different values of vscale
%   and epslevel.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

out = all(size(f.coeffs) == size(g.coeffs)) ...
    && all(f.values(:) == g.values(:)) ...
    && all(f.coeffs(:) == g.coeffs(:));

end