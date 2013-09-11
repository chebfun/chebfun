function [f, r] = addBreaksAtRoots(f, tol)
%ADDBREAKSATROOTS   Add breaks at appropriate roots of a CHEBFUN
%   ADDBREAKSATROOTS(F) introduces breakpoints at certain roots in the interior
%   of the domain of a CHEBFUN F. In particular, breaks are introduced at each
%   of the roots returned by ROOTS(F, 'nozerofun', 'nojump', 'noimps'), except
%   those which are deemed too close together or too close to existing
%   breakpoints.
%
%   ADDBREAKSATROOTS(F, TOL) provides a lower bound for the tolerance used in
%   the above exceptions.
%
%   If F is array-valued, breaks will be introduced in each of the columns at
%   unique(ROOTS(F)).
%
% See also ROOTS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Lower bound for tolerance:
if ( nargin == 1 )
    tol = 0;
end

% Locate roots:
r = roots(f, 'nozerofun', 'nojump', 'noimps');

% Since each column of an array-valued CHEBFUN must have the same breakpoints,
% we simply take unique(r(:)) and remove any remaining NaNs.
r = unique(r(:));
r(isnan(r)) = [];

% Discard any roots which are closer than the accuracy of the CHEBFUN:
el = epslevel(f);
hs = hscale(f);
rtol1 = max(el*hs, tol);
r([false ; diff(r) < rtol1]) = [];

% Avoid introducing new breakpoints close to an existing ones:
rtol2 = max(100*el*max(min(diff(f.domain)), 1), tol);
r(any(abs(bsxfun(@minus, r, f.domain)) < rtol2, 2)) = [];

% Add new breaks if required:
if ( ~isempty(r) )
    % Get the domain with the new breakpoints: (union is not required, by above)
    dom = unique([f.domain, r.']);

    % Introduce these breakpoints into f:
    f = restrict(f, dom);
end

end
