function A = cleanRows(A)
%CLEANROWS   Remove trailing zero rows from a matrix.
%   A = CLEANROWS(A) removes rows at the bottom of A which have all
%   entries negligible. A is typically the matrix of impulses.
%
% See also MERGECOLUMNS, MERGEIMPULSES, CLEANROWS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.
% Remove trivial columns:

if ( isempty(A) )
    return
end

% Get the tolerance:
tol = deltafun.pref.deltafun.deltaTol;

while( max(abs(A(end, :))) < tol )
    A(end, :) = [];
end