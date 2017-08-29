function W = fevalt(f, x, y, z)
%FEVALT   Evaluate a CHEBFUN3.
%   W = FEVALT(F, X, Y, Z) returns a tensor of values W of size 
%   length(X)-by-length(Y)-by-length(Z). X, Y and Z should be vectors of 
%   doubles. This is equivalent to fist constructing
%   [X,Y,Z] = NDGRID(X,Y,Z) and then computing W = FEVAL(F, X, Y, Z).
%   
%   FEVALT is different from W = FEVAL(F, X, Y, Z) which creates a
%   vector W while W = FEVALT(F, X, Y, Z) returns W as a tensor.
%
% See also CHEBFUN3/FEVAL.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    W = {[]}; 
    return
end

% Get low rank representation of f. 
[fCore,fCols, fRows, fTubes] = tucker(f);

colVals = feval(fCols, x(:));
rowVals = feval(fRows, y(:));
tubeVals = feval(fTubes, z(:));
W = chebfun3.txm(chebfun3.txm(chebfun3.txm(fCore, colVals, 1), ...
    rowVals, 2), tubeVals, 3);
end