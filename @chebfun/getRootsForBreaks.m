function [rBreaks, rAll] = getRootsForBreaks(f, tol)
%GETROOTSFORBREAKS   Get roots of a CHEBFUN and polish for use as breakpoints.
%   GETROOTSFORBREKAS(F) computes the roots of F (more specifically, ROOTS(F,
%   'nozerofun', 'nojump', 'nobreaks')) and then eliminates ones which are too
%   close together to be introduced into F (or into some other CHEBFUN with the
%   same domain as F) as breakpoints.  The roots are returned in a form
%   suitable for passing to ADDBREAKS().
%
%   GETROOTSFORBREAKS(F, TOL) provides a lower bound for the tolerance used in
%   deciding if two roots are too close.
%
%   [RBREAKS, RALL] = GETROOTSFORBREAKS(...) returns both the roots RBREAKS
%   deemed suitable for use as breakpoints and all of the roots computed by the
%   call to ROOTS.
%
% See also ROOTS, ADDBREAKS, ADDBREAKSATROOTS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Lower bound for tolerance:
if ( nargin == 1 )
    tol = 0;
end

% Locate roots:
rAll = roots(f, 'nozerofun', 'nojump', 'nobreaks');

% Reshape to a column vector, sort, and remove NaNs:
rBreaks = sort(rAll(:));
rBreaks(isnan(rBreaks)) = [];

% Discard any roots which are closer than the accuracy of the CHEBFUN (NB:
% This requires the roots to be sorted first.):
rootTol = max(eps*hscale(f), tol);
rBreaks([false ; diff(rBreaks) < rootTol]) = [];

end
