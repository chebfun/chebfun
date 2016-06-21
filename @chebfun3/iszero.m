function out = iszero(f)
%ISZERO   Check if a CHEBFUN3 is identically zero on its domain.
%   OUT = ISZERO(F) returns logical 1 if the CHEBFUN3 object F is exactly 
%   the zero function, and logical 0 otherwise.
%
% See also CHEBFUN2/ISZERO.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get data:
[fCore, fCols, fRows, fTubes] = tucker(f);

% Trivial check: If all the pivots are zero, then the SEPARABLEAPPROX is zero: 
if ( norm(fCore(:), inf) == 0 ) 
    out = 1; 
    return 
end

% Quick check: Evaluate on a grid. If the tensor is nonzero then the
% CHEBFUN3 is nonzero.
dom = f.domain; 
x = linspace(dom(1), dom(2), 10); 
y = linspace(dom(3), dom(4), 10);
z = linspace(dom(5), dom(6), 10);
vals = fevalt(f, x, y, z); 

if ( norm(vals(:), inf) > 0 ) 
   out = 0;
   return
end

% Slower check: The core may be nonzero, but the columns, rows, or tubes 
% may be zero:
[rk1, rk2, rk3] = rank(f);
out_cols = zeros(rk1, 1);
out_rows = zeros(rk2, 1);
out_tubes = zeros(rk3, 1);
for j = 1:rk1
    out_cols(j) = iszero(fCols(:, j));
end
for j = 1:rk2
    out_rows(j) = iszero(fRows(:, j));
end
for j = 1:rk3
    out_tubes(j) = iszero(fTubes(:, j));
end
out = ( all(out_cols) || all(out_rows) || all(out_tubes) );

end