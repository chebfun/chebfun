function z = fevalm(f, x, y)
% FEVALM   Evaluate a chebfun2

if ( isempty(f) )
    varargout = {[]}; 
    return
end

% Get the low rank representation for f. 
cols = f.cols; 
rows = f.rows; 
piv = f.pivotValues; 
d = 1./piv; 
d(d==inf) = 0;  % set infinite values to zero. 
zCol = feval(f.cols, y(:));
zRow = feval(f.rows, x(:));

z = zCol*diag( d )*zRow.';

end