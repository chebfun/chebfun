function a = all(f, pref)
%ALL    True if all elements of a chebfun are a nonzero number

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( nargin < 2 )
    pref = chebfun.pref();
end

% Check the impulses first (as this is trivial):
vs = get(f, 'vscale');
vs = max([vs{:}]);
el = get(f, 'epslevel');
tol = max(100*el*vs*pref.chebfun.eps/eps, eps);
a = all( f.impulses(1,:) > tol );

% Check to see if there are any roots:
if ( a == true )
    a = ~any(roots(f, 'nojump'));
end

end