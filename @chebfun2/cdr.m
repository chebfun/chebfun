function varargout = cdr( f )
%CDR decomposition of a chebfun2.
%
% [C,D,R] = CDR(F) produces a diagonal matrix D of size length(F) by
% length(F) and quasimatrices C and R of size inf by length(F) such 
% that f(x,y) = C(y,:) * D * R(x,:)'.
%
% D = CDR(F) returns a vector containing the pivot values used in the
% construction of F. 
%
% See also PIVOTS, SVD. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( isempty( f ) )
    varargout = {[]}; 
    return
end

% Get the low rank representation for f. 
cols = f.cols; 
rows = f.rows; 
piv = f.pivotValues; 
d = 1./piv; 
d( d == inf ) = 0;  % set infinite values to zero. 

% Output the CDR decomposition
if ( nargout <= 1 )
    varargout = { d };
else
    varargout = {cols, diag(d), rows};  % CDR decomposition
end

end