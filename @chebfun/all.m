function a = all(f)
%ALL   True if all elements of a CHEBFUN are a nonzero number.
%
% See also ANY, ISZERO.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check the impulses first (as this is trivial):
tol = vscale(f)*epslevel(f);
a = all(abs(f.impulses(:,:,1)) > tol, 1);

% Check to see if there are any roots:
if ( any(a == true) )
    r = roots(f, 'nojump');
    a = a & ~any(r, 1);
end

end
