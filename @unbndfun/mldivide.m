function X = mldivide(A, B)
%\   Left matrix divide for UNBNDFUN objects.
%   A\B returns the least-squares solution (with respect to the continuous L^2
%   norm) of A*X = B where A and B are UNBNDFUN objects.
%
%   The UNBNDFUN objects A and B are assumed to have the same domain. The method
%   gives no warning if their domains don't agree, but the output of the method
%   will be meaningless.
%
% See also QR, MRDIVIDE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call MLDIVIDE of the onefun fields of A and B.
X = A.onefun\B.onefun;

% % Alternatively we could call QR() at the UNBNDFUN level, since 
% %  A*X = (Q*R)*X = B ==> R*X = Q'*B ==> X = R\(Q'*B) = R\innerProduct(Q, B)
% which gives
% % Compute QR factorisation of A:
% % [Q, R] = qr(A, 0);
% 
% % Compute X:
% % X = R\innerProduct(Q, B);

end
