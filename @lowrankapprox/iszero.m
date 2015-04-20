function out = iszero( f )
%ISZERO   Check if a CHEBFUN2 is identically zero on its domain.
%   OUT = ISZERO( F ) return 1 if the CHEBFUN2 is exactly the zero function, and
%   0 otherwise. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get data:
pivots = f.pivotValues;
cols = f.cols;
rows = f.rows;

% Trivial check: If all the pivots are zero, then the CHEBFUN2 is zero: 
if ( norm(pivots, inf) == 0 ) 
    out = 1; 
    return 
end

% Quick check: Evaluate on a meshgrid. If the matrix is nonzero then the
% CHEBFUN2 is nonzero.
dom = f.domain; 
x = linspace(dom(1), dom(2), 10); 
y = linspace(dom(3), dom(4), 10);
[xx, yy] = meshgrid(x, y);
vals = feval(f, xx, yy); 
if ( norm( vals, inf ) > 0 ) 
   out = 0; 
   return
end

% Slower check: A pivot may be positive, but the columns or rows may be zero:
rk = length( pivots );
out_cols = zeros( length( pivots ), 1 );
out_rows = zeros( length( pivots ), 1 );
for j = 1:rk
    out_cols( j ) = iszero( cols(:, j) );
    out_rows( j ) = iszero( rows(:, j) );
end
bolslices = ( all(out_cols) || all(out_rows) );
out = bolslices;

end
