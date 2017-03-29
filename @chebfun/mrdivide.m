function X = mrdivide(B, A, flag)
%/   Right matrix divide for CHEBFUN objects.
%   X = B/A or X = mrdivide(B, A) is equivalent to X = B./A.
%
%   X = mrdivide(B, A, 'least-squares') or mrdivide(B, A, 'ls')returns the
%   least squares solution to X*A = B.
%
% See also RDIVIDE, MLDIVIDE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Parse inputs:
if ( nargin < 3 )
    flag = 0;
else
    flag = any(strcmpi(flag, {'least-squares', 'ls'}));
end

if ( ~flag ) % B./A
    X = rdivide(B, A);
else         % Least-squares solution
    X = mrdivide_ls(B, A);    
end

end

function X = mrdivide_ls(B, A)
% TODO: Give an example for each of the cases below.

if ( isscalar(A) )
    if ( A == 0 )
        % TODO:  Return identically Inf/NaN CHEBFUN instead?
        error('CHEBFUN:CHEBFUN:mrdivide:divisionByZero', ...
            'Division by zero.')
    end
    X = B*(1/A);

elseif ( size(A, 2) ~= size(B, 2) )
    error('CHEBFUN:CHEBFUN:mrdivide:dimensions', ...
        'Matrix dimensions must agree.')
    
elseif ( isnumeric(A) )
    % [INF x M] * [M x N] = [INF x N]:
    [Q, R] = qr(B, 0);
    X = Q * mrdivide(R,A);
    % X = B*(eye(size(B,2))/A);
    
elseif ( ~A(1).isTransposed )
    % [M x INF] * [INF x N] = [M x N]:
    [Q, R] = qr(A, 0);
    X = mrdivide_ls(B, R) * Q';
   
else
    % [M x N] * [N x INF] = [M x INF]:
    % XA = B
    % A^* X^* = B^*
    % QR X^*  = B^*
    % X R^* Q^* = B
    % X = (BQ)R^{-*}
    [Q, R] = qr(A', 0);
    X = mrdivide_ls((B*Q), R');

end

end
