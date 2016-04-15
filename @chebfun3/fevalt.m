function W = fevalt(f, x, y, z)
%FEVALT   Evaluate a CHEBFUN3.
% 
%   W = FEVALT(F, X, Y, Z) returns a tensor of values W of size 
%   length(X)-by-length(Y)-by-length(Z). X, Y, and Z should be vectors of 
%   doubles. This is equivalent to 
%   [x,y,z] = ndgrid(x,y,z) and then W = feval(f, x, y, z).

if ( isempty(f) )
    W = {[]}; 
    return
end

% Get the low rank representation for f. 
[fCore,fCols, fRows, fTubes] = st(f);

colVals = feval(fCols, x(:));
rowVals = feval(fRows, y(:));
tubeVals = feval(fTubes, z(:));
W = chebfun3.txm(chebfun3.txm(chebfun3.txm(fCore, colVals, 1), ...
    rowVals, 2), tubeVals, 3);
end