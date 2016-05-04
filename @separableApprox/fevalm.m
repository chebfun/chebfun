function out = fevalm(f, x, y)
% FEVALM   Evaluate a SEPARABLEAPPROX.
% 
% Z = FEVALM(F, X, Y) returns a matrix of values Z of size length(X)-by-length(Y). 
% X and Y should be vectors of doubles. This is equivalent to making a meshgrid 
% of the vectors X and Y and then using FEVAL to evaluate at that grid.
% 
% See also separableApprox/feval.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check for empty separableApprox object: 
if ( isempty( f ) )
    out = []; 
    return
end

% Get the low rank representation for f. 
[cols, D, rows] = cdr( f );

% Evaluate individual columns: 
zCol = feval( cols, y(:) );
zRow = feval( rows, x(:) );

% Equivalent to 
% [xx, yy] = meshgrid( x(:), y(:))
% out = feval(f, xx, yy)
out = zCol * D * zRow.';

end
