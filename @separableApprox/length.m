function [out1, out2] = length( F )
%LENGTH  The rank of a SEPARABLEAPPROX.
%   K = LENGTH(F) returns the rank of the SEPARABLEAPPROX representation.
%
%   [M, N] = LENGTH( F ) returns the length of the column and 
%   row slices employed in the separable model.  The interpretation 
%   of this quantity depends on the underlying representation. 
%
% See also RANK.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( F ) ) 
    out1 = []; 
    out2 = []; 
    return
end

if ( iszero( F ) ) 
    out1 = 1; 
    out2 = 1;
    return
end

% Extract length of underlying objects. 
if ( nargout <= 1 )
    out1 = length( F.pivotValues );
else
    out1 = length( F.rows );
    out2 = length( F.cols );
end

end
