function [out1, out2] = length(F)
%LENGTH  The rank of a chebfun2.
%
% K = LENGTH(F) returns the rank of the chebfun2.
%
% [m, n] = LENGTH(F) returns the polynomial degree of the column and row
%        slices.
%
% See also RANK.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( nargout <= 1 )
    out1 = length( F.pivotValues );
else
    out1 = length( F.rows );
    out2 = length( F.cols );
end

end