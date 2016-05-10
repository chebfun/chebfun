function varargout = cdr( f )
%CDR decomposition of a SEPARABLEAPPROX.
%   [C,D,R] = CDR(F) produces a diagonal matrix D of size length(F) by length(F)
%   and quasimatrices C and R of size inf by length(F) such that f(x,y) = C(y,:)
%   * D * R(x,:)'.
%
%   D = CDR(F) returns a vector containing the pivot values used in the
%   construction of F.
%
% See also PIVOTS, SVD. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty( f ) )
    varargout = cell(1, nargout); 
    return
end

% Get the low rank representation for f. 
cols = f.cols; 
rows = f.rows; 
piv = f.pivotValues; 
d = 1./piv; 

% Set infinite values to zero.
d( d == inf ) = 0;  

% Output:
if ( nargout <= 1 )
    varargout = { d };
else
    % CDR decomposition
    varargout = {cols, diag(d), rows};  
end

end
