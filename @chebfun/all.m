function a = all(f)
%ALL   True if all values of a CHEBFUN are nonzero.
%
% See also ANY, ISZERO.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check the piont values first (as this is trivial):
tol = vscale(f)*eps;
a = all(abs(f.pointValues) > tol);

% Check to see if there are any roots:
if ( any(a) )
    r = roots(f, 'nojump');
    if ( isempty(r) )
        return
    end
    
    indNaN = isnan(r);
    noRoots = all(indNaN, 1);
    a = a & noRoots;
end

end
