function X = mldivide(A, B)
%\   Left matrix divide for FOURTECH objects.
%   A\B returns the least-squares solution (with respect to the continuous L^2
%   norm) of A*X = B where A and B are FOURTECH objects.
%
% See also QR, MRDIVIDE.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ~isa(A, 'fourtech') || ~isa(B, 'fourtech') || ~strcmp(class(A), class(B)) )
    error('CHEBFUN:FOURTECH:mldivide:fourtechMldivideUnknown', ...
            [ 'Arguments to FOURTECH mldivide must both be FOURTECH objects' ...
              'of the same type.' ]);
end

% NOTE: A*X = (Q*R)*X = B ==> R*X = Q'*B ==> X = R\(Q'*B) = R\innerProduct(Q, B)

% Compute QR factorisation of A:
[Q, R] = qr(A, 0);

% Compute X:
X = R\innerProduct(Q, B);

end