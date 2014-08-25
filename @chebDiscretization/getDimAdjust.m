function dimAdjust = getDimAdjust(L)
%GETDIMADJUST   Adjust dimension of discretization.
%   GETDIMADJUST(L) returns zero unless L is a LINOP.
%
%   A = GETDIMADJUST(L), where L is a LINOP,  returns a matrix of values, A, of
%   the same dimension as L.BLOCKS, which informs discretization methods how
%   the dimension L.DIM of that block must be adjusted (usually increased) so
%   that the highest order derivatives of each of the variables appearing in
%   the LINOP system are discretized at the same dimension.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ~isa(L, 'linop') )
    
    % No adjustment
    dimAdjust = 0;
    
else
    
    % The input adjustment size of the (j,k) entry is max(diffOrder(:,k))
    dimAdjust = max(getDiffOrder(L), [], 1);
    dimAdjust = max(dimAdjust, 0);
    dimAdjust = repmat(dimAdjust, size(L, 1), 1);
    
end
    
end
