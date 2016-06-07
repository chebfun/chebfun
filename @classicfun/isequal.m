function out = isequal(f, g)
%ISEQUAL   Test if two CLASSICFUN objects are equal.
%   ISEQUAL(F, G) returns logical 1 (TRUE) if the CLASSICFUN objects F and G have the
%   same domain and .onefun fields which are equal, as determined by their
%   respective class methods.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check whether .domains and .onefun fields agree:
out = (all(f.domain == g.domain) && isequal(f.onefun, g.onefun));

end
