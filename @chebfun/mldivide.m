function X = mldivide(A, B)
%\   Left matrix divide for CHEBFUN objects.
%   A\B in general gives the least squares solution to A*X = B.
%
% See also MRDIVIDE, LDIVIDE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Something is wrong with the comment lines below.
% TODO: Give an example for each of the cases below.

if ( isscalar(A) )
    % Trivial case (division by a constant):
    X = (A\1) * B;
    
elseif ( size(A, 1) ~= size(B, 1) )
    error('CHEBFUN:CHEBFUN:mldivide:agree', 'Matrix dimensions must agree.')
        
elseif ( isnumeric(A) )
    % [M x N] * [N x INF] = [M x INF]:
    
    [Q, R] = qr(B', 0);
    X = (A\R') * Q';
    % X = (A\eye(size(B,1)))*B;
    
elseif ( A(1).isTransposed )
    % [M x INF] * [INF x N] = [M x N]:
    %        AX = B
    %   X^* A^* = B^*
    %   X^* QR  = B^*
    % R^* Q^* X = B
    %         X = Q(R^{-*} B)
    
    [Q, R] = qr(A', 0);
    X = Q * (R'\B);
    
else
    % [INF x N] * [N x M] = [INF x M]:
    
    [Q, R] = qr(A, 0);
    X = R \ innerProduct(Q, B);
    
end
