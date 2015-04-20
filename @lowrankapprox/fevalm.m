function z = fevalm(f, x, y)
% FEVALM   Evaluate a CHEBFUN2.
% 
% Z = FEVALM(F, X, Y) returns a matrix of values Z of size length(X)-by-length(Y). 
% X and Y should be vectors of doubles. This is equivalent to making a meshgrid 
% of the vectors X and Y and then using FEVAL to evaluate at that grid.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    varargout = {[]}; 
    return
end

% Get the low rank representation for f. 
[cols, D, rows] = cdr(f);

zCol = feval(f.cols, y(:));
zRow = feval(f.rows, x(:));
z = zCol * D * zRow.';

end
