function A = cleanRows(A, pref)
%CLEANROWS   Remove trailing zero rows from a matrix.
%   A = CLEANROWS(A) removes rows at the bottom of A which have negligible
%   entries. A is typically the matrix of impulses. This function uses the
%   tolerance provided by CHEBFUNPREF and uses that tolerance to decide whether
%   a delta function is trivial or not.
%
% See also MERGECOLUMNS, MERGEDELTAS, CLEANCOLUMNS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(A) )
    return
end

% Get the tolerance:
if ( nargin < 3 || isempty(pref) )
    pref = chebfunpref();
end
deltaTol = pref.deltaPrefs.deltaTol;

% Remove trivial rows:
while ( max(abs(A(end,:))) < deltaTol )
    A(end,:) = [];
end

end
