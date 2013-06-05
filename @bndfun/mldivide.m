function X = mldivide(A, B)
%\   Left matrix divide for BNDFUN objects.
%
%   A\B returns the least-squares solution (with respect to the continuous L^2
%   norm) of A*X = B where A and B are BNDFUN objects.
%
% See also QR, MRDIVIDE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% NOTE: A*X = (Q*R)*X = B ==> R*X = Q'*B ==> X = R\(Q'*B) = R\innerProduct(Q, B)
% Compute QR factorisation of A:
[Q, R] = qr(A, 0);

% Compute X:
X = R\innerProduct(Q, B);

% TODO: Do we want to compute MLDIVIDE using the method above, or call MLDIVIDE
% at the ONEFUN level, e.g. via
% X = A.onefun\B.onefun;


end
