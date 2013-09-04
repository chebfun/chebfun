function f = addBreaksAtRoots(f)
% Insert breaks whre the function f is zero (with some tolerance).

% [TODO]: This is used in other places (such as ABS) and needs de-duplicating.

% Locate roots:
r = roots(f, 'nozerofun', 'nojump', 'noimps');

% Since each column of an array-valued CHEBFUN must have the same breakpoints,
% we simply take unique(r(:)) and remove any remaining NaNs.
r = unique(r(:));
r(isnan(r)) = [];

% Discard any roots which are closer than the accuracy of the CHEBFUN:
el = epslevel(f);
hs = hscale(f);
rtol1 = el*hs;
r([false ; diff(r) < rtol1]) = [];

% Avoid introducing new breakpoints close to an existing ones:
rtol2 = 100*el*max(min(diff(f.domain)), 1);
r(any(abs(bsxfun(@minus, r, f.domain)) < rtol2, 2)) = [];

% Add new breaks if required:
if ( ~isempty(r) )
    % Get the domain with the new breakpoints: (union is not required, by above)
    dom = unique([f.domain, r.']);

    % Introduce these breakpoints to f:
    f = restrict(f, dom);
end

end

