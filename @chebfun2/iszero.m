function out = iszero( f )
% ISZERO   Check if a CHEBFUN2 is identically zero on its domain.
%
% OUT = ISZERO( F ) return 1 if the CHEBFUN2 is exactly the zero function, and
% 0 otherwise. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

pivots = f.pivotValues;
cols = f.cols;
rows = f.rows;

rk = length( pivots );
out_pivots = ( norm(pivots, inf) == 0 );
out_cols = zeros( length( pivots ), 1 );
out_rows = zeros( length( pivots ), 1 );
for j = 1:rk
    out_cols( j ) = iszero( cols(:, j) );
    out_rows( j ) = iszero( rows(:, j) );
end

bolslices = ( all(out_cols) || all(out_rows) );

out = bolslices || out_pivots;

end
