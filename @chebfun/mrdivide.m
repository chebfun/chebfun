function C = mrdivide(B, A)
%/   Right matrix divide for CHEBFUN objects.
%   B/A in general gives the least squares solution to X*A = B.
%
% See also RDIVIDE, MLDIVIDE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( isscalar(A) )
    C = B*(1/A);
    
elseif ( size(A, 2) ~= size(B, 2) )
    error('CHEBFUN:mldivide:dimensions', 'Matrix dimensions must agree.')
    
elseif ( isnumeric(A) )
    % C = B*(eye(size(B,2))/A);
    [Q, R] = qr(B, 0);
    C = Q * (R/A);
    
elseif ( ~A.isTransposed )
    [Q,R] = qr(A, 0);
    C = (B/R) * Q';
    
else
    [Q,R] = qr(A', 0);
    C = (B*Q) / R';
    
end
