function [A, v] = cleanColumns(A, v)
%CLEANCOLUMNS   Remove zero columns from a matrix.
%   [A v] = CLEANCOLUMNS(A, V) removes columns which have all entries
%   negligible and removes the corresponding entry in the vector V. A is
%   typically the matrix of impulses and V is the vector containing
%   locations.
%
% See also MERGECOLUMNS, MERGEIMPULSES, CLEANROWS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.
% Remove trivial columns:

tol = deltafun.pref.deltafun.deltaTol;

m = size(A, 2);
if ( length(v) ~= m || size(v, 1) > 1 )
    error( 'CHEBFUN:DELTAFUN:cleanColumns', 'No. of columns of A should equal the length of the vector va' );
end

j = 1;
for k = 1:m
    if( max(abs(A(:, j))) <  tol )
        v(j) = [];
        A(:, j) = [];
    else
        j = j + 1;
    end
end