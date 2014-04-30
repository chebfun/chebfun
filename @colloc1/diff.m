function D = diff(disc, m)
%DIFF    Differentiation operator for COLLOC1 discretization.
%   D = DIFF(DISC) gives the matrix such that if v=D*u, then v=u', where u
%   is a COLLOC1 representation of a Chebyshev polynomial.
%
%   DIFF(DISC, M) for positive integer M returns D^M (through a better
%   algorithm than multiplication).

%  Copyright 2014 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% Store information about domain and dimensions.
d = disc.domain;
n = disc.dimension;

if m == 0
    % Trivial case
    D = eye(sum(n));
else
    numIntervals = disc.numIntervals;
    
    % Find the diagonal blocks.
    blocks = cell(numIntervals);
    for k = 1:numIntervals
        len = d(k+1) - d(k);
        blocks{k} = diffmat(n(k),m) * (2/len)^m; % Scaled diffmats
    end
    
    % Assemble!
    D = blkdiag(blocks{:});
end

end

function D = diffmat(n, k)
%DIFFMAT   Chebyshev differentiation matrix
%   D = DIFFMAT(N) is the matrix that maps function values at N Chebyshev
%   points to values of the derivative of the interpolating polynomial at
%   those points.
%
%   D = DIFFMAT(N,K) is the same, but for the Kth derivative.

% TODO: Duplicated?
% TODO: Cache this?

if ( nargin < 2 ), k = 1; end
if ( n == 0 ), D = []; return, end
if ( n == 1 ), D = 0; return, end

% 1st-kind Chebyshev grid:
x = chebtech1.chebpts(n);
% 1st-kind Barycentric weights:
w = chebtech1.barywts(n);

ii = (1:n+1:n^2)';              % indices of diagonal
Dx = bsxfun(@minus,x,x');       % all pairwise differences
Dx(ii) = Dx(ii) + 1;            % add identity
Dxi = 1./Dx;                    % reciprocal 
Dw = bsxfun(@rdivide,w.',w);    % pairwise divisions
Dw(ii) = Dw(ii) - 1;            % subtract identity

% k = 1
D = Dw .* Dxi;
D(ii) = 0; D(ii) = - sum(D,2);                  % negative sum trick

if ( k > 1 )
    % 2nd order
    D = 2*D .* (repmat(D(ii),1,n) - Dxi);
    D(ii) = 0; D(ii) = - sum(D,2);              % negative sum trick

    % higher orders
    for m = 3:k
        D = m*Dxi .* (Dw.*repmat(D(ii),1,N) - D);
        D(ii) = 0; D(ii) = - sum(D,2);          % negative sum trick
    end
end

end

