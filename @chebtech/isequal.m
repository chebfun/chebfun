function out = isequal(f, g)
%ISEQUAL Test if CHEBTECH objects are equal.
%    ISEQUAL(F, G) returns logical 1 (TRUE) if the CHEBTECH objects F and G are
%    the same length and contain the same values. They may have different 
%    values of vscale and epslevel.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

out = all(size(f.values) == size(g.values)) ... % Size of values is the same
    && all(f.values(:) == g.values(:)) ...      % Entries of values are same
    && all(f.coeffs(:) == g.coeffs(:));         % Coefficients are the same.

end
