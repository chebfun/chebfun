function X = mrdivide(B, A)
%/   Right matrix divide for CHEBFUN objects.
%   B/A in general gives the least squares solution to X*A = B.
%
% See also RDIVIDE, MLDIVIDE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( isscalar(A) )
    X = B*(1/A);
    
elseif ( size(A, 2) ~= size(B, 2) )
    error('CHEBFUN:mldivide:dimensions', 'Matrix dimensions must agree.')
    
elseif ( isnumeric(A) )
    % [INF x M] * [M x N] = [INF x N]:
    [Q, R] = qr(B, 0);
    X = Q * (R/A);
    % X = B*(eye(size(B,2))/A);
    
elseif ( ~A.isTransposed )
    % [M x INF] * [INF x N] = [M x N]:
    [Q, R] = qr(A, 0);
    X = (B/R) * Q'; % TODO: .'?
    
else
    % [M x N] * [N x INF] = [M x INF]:
    [Q, R] = qr(A', 0); % TODO: .'?
    X = (B*Q) / R';

end
