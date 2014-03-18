function X = mldivide(A, B)
%\   Left matrix divide for FOURIERTECH objects.
%
%   A\B returns the least-squares solution (with respect to the continuous L^2
%   norm) of A*X = B where A and B are FOURIERTECH objects.
%
% See also QR, MRDIVIDE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ~isa(A, 'fouriertech') || ~isa(B, 'fouriertech') || ~strcmp(class(A), class(B)) )
    error('CHEBFUN:FOURIERTECH:mldivide:chebtechMldivideUnknown', ...
            [ 'Arguments to FOURIERTECH mldivide must both be FOURIERTECH objects' ...
              'of the same type.' ]);
end

error('CHEBFUN:FOURIERTECH:mldivide','This option requires the QR which is not implemented for FOURIERTECH objects yet');

% % NOTE: A*X = (Q*R)*X = B ==> R*X = Q'*B ==> X = R\(Q'*B) = R\innerProduct(Q, B)
% 
% % Compute QR factorisation of A:
% [Q, R] = qr(A, 0);
% 
% % Compute X:
% X = R\innerProduct(Q, B);

end
