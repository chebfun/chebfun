function out = isequal(f, g)
%ISEQUAL   Test if BNDFUN objects are equal.
%   ISEQUAL(F, G) returns logical 1 (TRUE) if the BNDFUN objects F and G have
%   onefun fields which are equal, as determined by the respective class
%   methods.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% If their domains don't agree, they're surely not equal!
if ( ~checkDomain(f, g))
    out = 0;
else
    % Check whether the onefun fields of f and g are equal
    out = isequal(f.onefun, g.onefun);
end

end
