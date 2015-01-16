function [out1, out2] = length(F)
%LENGTH  The rank of a CHEBFUN2.
%   K = LENGTH(F) returns the rank of the CHEBFUN2.
%
%   [M, N] = LENGTH(F) returns the polynomial degree of the column and row
%   slices.
%
% See also RANK.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(F) ) 
    out1 = []; 
    out2 = []; 
    return
end

if ( iszero(F) ) 
    out1 = 1; 
    out2 = 1;
    return
end

if ( nargout <= 1 )
    out1 = length(F.pivotValues);
else
    out1 = length(F.rows);
    out2 = length(F.cols);
end

end
