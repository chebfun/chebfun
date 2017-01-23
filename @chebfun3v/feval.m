function vals = feval(F, x, y, z)
%FEVAL   Pointwise evaluation of a CHEBFUN3V object.
%   F(X, Y, Z) returns value of F at the point (X,Y,Z).
%
% See also CHEBFUN3V/SUBSREF.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(F) )
    vals = []; 
    return
end

nF = F.nComponents;
if ( isnumeric(x) && isnumeric(y) && isnumeric(z) )
    
    vals = zeros(nF, length(x));
    % Evaluate each component:
    for jj = 1:nF
        vals(jj, :) = feval(F.components{jj}, x, y, z);  
    end

elseif ( isnumeric(x) && strcmpi(y, ':') && strcmpi(z, ':') ||...
        strcmpi(x, ':') && isnumeric(y) && strcmpi(z, ':') ||...
        strcmpi(x, ':') && strcmpi(y, ':') && isnumeric(z) )
    % Is the output supposed to be a chebfun2v object?
    
    vals = chebfun2v();
    vals.nComponents = nF;
    vals.isTransposed = F.isTransposed;
    % Evaluate each component:
    for jj = 1:nF
        vals.components{jj} = feval(F.components{jj}, x, y, z);
    end
    
elseif( isnumeric(x) && isnumeric(y) &&  strcmpi(z, ':') || ...
        isnumeric(x) && strcmpi(y, ':') && isnumeric(z) || ...
        strcmpi(x, ':') && isnumeric(y) && isnumeric(z) )
    % Is the output supposed to be a vector of 1D chebfuns, i.e., a 
    % quasimatrix?
    vals = chebfun(zeros(length(x), nF));
    for jj = 1:nF
        vals(:, jj) = feval(F.components{jj}, x, y, z);
    end
    
end

end