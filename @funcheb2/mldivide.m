function X = mldivide(A, B)
%\	Left matrix divide for FUNCHEB2 objects.
%
%   A\B returns the least squares solution (wrt the continuous L^2 norm) of A*X
%   = B where A and B FUNCHEB2 objects.
%
% See also QR, MRDIVIDE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ~isa(A, 'funcheb2') || ~isa(B, 'funcheb2') )
    error('CHEBFUN:FUNCHEB2:mldivide:funcheb2MldivideUnknown', ...
        'Arguments to FUNCHEB2 mldivide must both be FUNCHEB2 objects.');
end

% NOTE: A*X = (Q*R)*X = B ==> R*X = Q'*B ==> X = R*(Q'*B) = R*innerProduct(Q, B)

% Compute QR factorisation of A:
[Q, R] = qr(A, 0);

% Compute x:
X = R\innerProduct(Q, B);

end
