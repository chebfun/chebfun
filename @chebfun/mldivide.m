function C = mldivide(A, B)
%\   Left matrix divide.
%   A\B in general gives the least squares solution to A*X = B.
%
% See also MRDIVIDE, LDIVIDE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( isscalar(A) )
    C = (A\1) * B;
    
elseif ( size(A, 1) ~= size(B, 1) )
    error('CHEBFUN:mldivide:agree', 'Matrix dimensions must agree.')
    
elseif ( isnumeric(A) )
    % C = (A\eye(size(B,1)))*B;
    [Q,R] = qr(B',0);
    C = (A\R') * Q';
    
elseif ( A.isTransposed )
    [Q,R] = qr(A',0);
    C = Q * (R'\B);
    
else
    [Q,R] = qr(A,0);
    C = R \ (Q'*B);
    
end