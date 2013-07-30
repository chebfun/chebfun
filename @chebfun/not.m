function f = not(f, pref)
%~   Chebfun logical NOT.
%   NOT(F) returns a chebfun which evaluates to zero at all points where F is zero
%   and one otherwise.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Loop over the funs:
for k = 1:numel(f.funs)
    f.funs{k} = not(f.funs{k});
end

if ( nargin < 2 )
    pref = chebfun.pref();
end

% Impulses:
vs = get(f, 'vscale'); vs = max([vs{:}]);
tol = max(100*pref.chebfun.eps*vs, eps);
f.impulses = ~(f.impulses(1,:) < tol);

end
