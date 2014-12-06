function f = addBreaks(f, breaks, tol)
%ADDBREAKS   Add breakpoints to the domain of a CHEBFUN.
%   ADDBREAKS(F, BREAKS) attempts to insert breakpoints in F at the points in
%   the vector BREAKS.  BREAKS need not be sorted or have only unique values,
%   but it should consist only of points in the domain of F.  Breakpoints will
%   not be inserted if they are too close together or too close to existing
%   breakpoints.
%
%   ADDBREAKS(F, BREAKS, TOL) does the same but uses the tolerance TOL as a
%   lower bound for the tolerance used in deciding if breakpoints are too close
%   to each other or to existing ones to qualify for insertion.
%
% See also ADDBREAKSATROOTS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Lower bound for tolerance:
if ( nargin == 2 )
    tol = 0;
end

% Since each column of an array-valued CHEBFUN must have the same breakpoints,
% we simply take unique(breaks(:)).
breaks = unique(breaks(:));

% Get rid of new breakpoints at infinity:
breaks( isinf(breaks) ) = [];

% Get rid of new breakpoints at infinity:
domain = f.domain( isfinite(f.domain) );

% Set a tolerance:
breakTol = max(100*epslevel(f)*max(min(diff(domain)), 1), tol);

% Distance between new breakpoints and existing ones:
dist = abs(bsxfun(@minus, breaks, domain));

% Avoid introducing new breakpoints close to existing ones: 
% Note that: 1. In case of [-Inf Inf], DIST is empty.
%            2. In case of [a Inf] or [-Inf b], BREAKTOL is empty.
if ( ~isempty(dist) && ~isempty(breakTol) )
    breaks( any(dist < breakTol, 2) ) = [];
end

% Add new breaks if required:
if ( ~isempty(breaks) )
    % Get the domain with the new breakpoints: (union is not required, by above)
    dom = unique([f.domain, breaks.']);
    
    % Handle periodic functions:
    if ( isPeriodicTech(f) )
        f = chebfun(f);
    end

    % Introduce these breakpoints into f:
    f = restrict(f, dom);
end

end
