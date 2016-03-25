function f = simplify( f, tol ) 
% Simplify a SEPARABLEAPPROX
% 
% F = SIMPLIFY( F ) compressed the representation of F to one that is
% numerically the same, but requires fewer parameters to store. Currently this
% simplifies the polynomial degree of F, but not the rank. 
%
% F = SIMPLIFY(F, TOL) does the same as above but uses the scalar TOL instead
% of the default simplification tolerance as the relative threshold level.

% Copyright 2015 by The University of Oxford and The Chebfun2 Developers.
% See http://www.chebfun.org/ for Chebfun2 information.

if ( nargin < 2 )
    tol = []; 
end

% Simplify the column and row slices. 
f.cols = simplify( f.cols, tol, 'globaltol' ); 
f.rows = simplify( f.rows, tol, 'globaltol' ); 

% Ensure the left and right limits match the endpoints:
f.cols = resetPointValues(f.cols);
f.rows = resetPointValues(f.rows);

end 
