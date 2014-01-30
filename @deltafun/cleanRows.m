function A = cleanRows(A, pref)
%CLEANROWS   Remove trailing zero rows from a matrix.
%   A = CLEANROWS(A) removes rows at the bottom of A which have negligible
%   entries. A is typically the matrix of impulses.
%
% See also MERGECOLUMNS, MERGEIMPULSES, CLEANCOLUMNS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% TODO: Comment on tolerance.

if ( isempty(A) )
    return
end

% Get the tolerance:
if ( nargin < 3 || isempty(pref) )
    pref = chebpref();
end
deltaTol = pref.deltaPrefs.deltaTol;

% Remove trivial rows:
while ( max(abs(A(end,:))) < deltaTol )
    A(end,:) = [];
end

end