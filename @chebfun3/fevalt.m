function W = fevalt(f, x, y, z)
%FEVALT   Evaluate a CHEBFUN3.
% 
%   W = FEVALT(F, X, Y, Z) returns a tensor of values W of size 
%   length(X)-by-length(Y)-by-length(Z). X, Y, and Z should be vectors of 
%   doubles. This is equivalent to 
%   [x,y,z] = ndgrid(x,y,z) and then W = feval(f, x, y, z).
%   
%   The difference with w = feval(f, x, y, z) is that feval creates a
%   vector w, while w = feval(f, x, y, z) form w as a tensor.
%
%   See also chebfun3/feval.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    W = {[]}; 
    return
end

% Get the low rank representation for f. 
[fCore,fCols, fRows, fTubes] = tucker(f);

colVals = feval(fCols, x(:));
rowVals = feval(fRows, y(:));
tubeVals = feval(fTubes, z(:));
W = chebfun3.txm(chebfun3.txm(chebfun3.txm(fCore, colVals, 1), ...
    rowVals, 2), tubeVals, 3);
end