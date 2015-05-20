function f = addBreaksAtRoots(f, tol)
%ADDBREAKSATROOTS   Add breaks at appropriate roots of a CHEBFUN.
%   ADDBREAKSATROOTS(F) introduces breakpoints at certain roots in the interior
%   of the domain of a CHEBFUN F. In particular, breaks are introduced at each
%   of the roots returned by ROOTS(F, 'nozerofun', 'nojump', 'nobreaks'), except
%   those which are deemed too close together or too close to existing
%   breakpoints.
%
%   ADDBREAKSATROOTS(F, TOL) provides a lower bound for the tolerance used in
%   the above exceptions.
%
%   If F is array-valued, breaks will be introduced in each of the columns at
%   unique(ROOTS(F)).
%
% See also ADDBREAKS, ROOTS, GETROOTSFORBREAKS, DEFINEPOINT.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%TODO: return a quasimatrix from array-valued CHEBFUN input?

% Parse inputs:
if ( nargin == 1 )
    tol = 0;
elseif ( isa(tol, 'chebfunpref') )
    tol = tol.techPrefs.eps;
end

for k = 1:numel(f)
    f(k) = columnAddBreaksAtRoots(f(k), tol);
end

end

function f = columnAddBreaksAtRoots(f, tol)

% Get the roots:
[rBreaks, rAll] = getRootsForBreaks(f, tol);

% Add new breaks if required:
if ( ~isempty(rBreaks) )
    oldDomain = f.domain;
    f = addBreaks(f, rBreaks, tol);

    % Enforce zero impulses at roots only if new breakpoints were added (i.e.,
    % the roots were not too close to existing breakpoints):
    if ( ~isequal(f.domain, oldDomain) )
        for k = 1:min(size(f))
            % TODO: Allow a tolerance?
            f.pointValues(ismember(f.domain, rAll(:,k)), k, :) = 0;
        end
    end
end

end
