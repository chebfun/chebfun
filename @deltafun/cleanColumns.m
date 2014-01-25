function [A, v] = cleanColumns(A, v, pref)
%CLEANCOLUMNS   Remove zero columns from a matrix.
%   [A v] = CLEANCOLUMNS(A, V) removes columns which have all entries
%   negligible and removes the corresponding entry in the vector V. A is
%   typically the matrix of impulses and V is the vector containing
%   locations.
%
% See also MERGECOLUMNS, MERGEIMPULSES, CLEANROWS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Get the tolerance:
if ( nargin < 3 || isempty(pref) )
    pref = chebpref();
end
deltaTol = pref.deltaPrefs.deltaTol;

m = size(A, 2);
if ( length(v) ~= m || size(v, 1) > 1 )
    error( 'DELTAFUN:cleanColumns', 'No. of columns of A should equal the length of the vector v.' );
end

j = 1;
for k = 1:m
    if( max(abs(A(:, j))) <  deltaTol )
        v(j) = [];
        A(:, j) = [];
    else
        j = j + 1;
    end
end