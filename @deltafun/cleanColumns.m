function [A, v] = cleanColumns(A, v)
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