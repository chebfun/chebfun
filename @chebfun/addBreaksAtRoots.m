function f = addBreaksAtRoots(f, tol)
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
% See also ADDBREAKS, ROOTS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Lower bound for tolerance:
if ( nargin == 1 )
    tol = 0;
end

% Locate roots:
rAll = roots(f, 'nozerofun', 'nojump', 'noimps');

% Reshape to a column vector and remove NaNs:
r = rAll(:);
r(isnan(r)) = [];

% Discard any roots which are closer than the accuracy of the CHEBFUN:
rootTol = max(epslevel(f)*hscale(f), tol);
r([false ; diff(r) < rootTol]) = [];

% Add new breaks if required:
if ( ~isempty(r) )
    f = addBreaks(f, r, tol);

    % Enforce zero impulses at roots:
    for k = 1:min(size(f))
        % TODO: Allow a tolerance?
        f.impulses(ismember(f.domain, rAll(:,k)), k, :) = 0;
    end
end

end
