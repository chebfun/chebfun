function y = clenshaw(x, c)
%CLENSHAW   Clenshaw's algorithm for evaluating a Chebyshev polynomial.
%   If C is a column vector, Y = CLENSHAW(X, C) evaluates the Chebyshev
%   expansion
%
%     Y = P_N(X) = C(1)*T_N(X) + ... + C(N)*T_1(X) + C(N+1)*T_0(X)
%
%   using Clenshaw's algorithm.
%
%   If C is an (N+1) x M matrix, then CLENSHAW interprets each of the columns
%   of C as coefficients of a degree N polynomial and evaluates the M Chebyshev
%   expansions
%
%     Y_m = P_N(X) = C(1,m)*T_N(X) + ... + C(N,m)*T_1(X) + C(N+1,m)*T_0(X)
%
%   for 1 <= m <= M, returning the results as columns of a matrix Y =
%   [Y_1 ... Y_M].
%
%   In both cases, X must be a column vector.
%
% See also CHEBTECH.FEVAL, CHEBTECH.BARY.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer note: Clenshaw is not typically called directly, but by FEVAL().
%
% Developer note: The algorithm is implemented both for scalar and for vector
% inputs. Of course, the vector implementation could also be used for the scalar
% case, but the additional overheads make it a factor of 2-4 slower. Since the
% code is short, we live with the minor duplication.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% X should be a column vector.
if ( size(x, 2) > 1 )
    warning('CHEBFUN:CHEBTECH:clenshaw:xDim', ...
        'Evaluation points should be a column vector.');
    x = x(:);
end

% Scalar or array-valued function?
if ( size(c, 2) == 1 )
    y = clenshaw_scl(x, c);
else
    y = clenshaw_vec(x, c);
end

end

function y = clenshaw_scl(x, c)
% Clenshaw scheme for scalar-valued functions.
bk1 = zeros(size(x)); 
bk2 = bk1;
x = 2*x;
for k = 1:size(c,1)-1
    bk = c(k) + x.*bk1 - bk2;
    bk2 = bk1; 
    bk1 = bk;
end
y = c(end) + .5*x.*bk1 - bk2;
end

function y = clenshaw_vec(x, c)
% Clenshaw scheme for array-valued functions.
x = repmat(x(:), 1, size(c, 2));
bk1 = zeros(size(x, 1), size(c, 2)); 
bk2 = bk1;
e = ones(size(x, 1), 1);
x = 2*x;
for k = 1:size(c, 1)-1
    bk = e*c(k,:) + x.*bk1 - bk2;
    bk2 = bk1; 
    bk1 = bk;
end
y = e*c(end,:) + .5*x.*bk1 - bk2;
end
