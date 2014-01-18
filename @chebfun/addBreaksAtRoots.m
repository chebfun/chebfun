function f = addBreaksAtRoots(f, tol)
%ADDBREAKSATROOTS   Add breaks at appropriate roots of a CHEBFUN
%
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
% See also ADDBREAKS, ROOTS, GETROOTSFORBREAKS, DEFINEPOINT

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Lower bound for tolerance:
if ( nargin == 1 )
    tol = 0;
end

% Get the roots:
[rBreaks, rAll] = getRootsForBreaks(f, tol);

% Add new breaks if required:
if ( ~isempty(rBreaks) )
    f = addBreaks(f, rBreaks, tol);

    % Enforce zero impulses at roots:
    for k = 1:min(size(f))
        % TODO: Allow a tolerance?
        f.impulses(ismember(f.domain, rAll(:,k)), k, :) = 0;
    end
end

end
