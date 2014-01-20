function bol = iszero( f )
% ISZERO   Check if a chebfun2 is identically zero on its domain.

pivots = f.pivotValues;
cols = f.cols;
rows = f.rows;

rk = length( pivots );
bolpivots = ( norm(pivots, inf) == 0 );
bolcols = zeros( length( pivots ), 1 );
bolrows = zeros( length( pivots ), 1 );
for j = 1 : rk
    bolcols( j ) = iszero( cols(:, j) );
    bolrows( j ) = iszero( rows(:, j) );
end

bolslices = ( any(bolcols) || any(bolrows) );

bol = bolslices || bolpivots;
end