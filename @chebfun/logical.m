function f = logical(f, pref)
%LOGICAL   Chebfun logical.
% LOGICAL(F) returns a chebfun which evaluates to one at all points where F is
% non-zero and zero otherwise.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Loop over the funs:
for k = 1:numel(f.funs)
    f.funs{k} = logical(f.funs{k});
end

if ( nargin < 2 )
    pref = chebfun.pref();
end

% Impulses:
vs = get(f, 'vscale'); vs = max([vs{:}]);
tol = max(100*pref.chebfun.eps*vs, eps);
f.impulses = (f.impulses(1,:) < tol);

end
