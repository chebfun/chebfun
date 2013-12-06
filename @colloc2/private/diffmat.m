function D = diffmat(N,k)

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
% DIFFMAT  Chebyshev differentiation matrix
% D = DIFFMAT(N) is the matrix that maps function values at N Chebyshev
% points to values of the derivative of the interpolating polynomial at
% those points.
%
% D = DIFFMAT(N,K) is the same, but for the Kth derivative.
%
% The matrices are computed using the 'hybrid' formula of Schneider &
% Werner [1] and Welfert [2] proposed by Tee [3].

% Copyright 2011 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% References:
%  [1] Schneider, C. and Werner, W., "Some new aspects of rational
%   interpolation", Math. Comp. (47) 285--299, 1986.
%  [2] Welfert, B. D., "Generation of pseudospectral matrices I", SINUM,
%   (34) 1640--1657.
%  [3] Tee, T. W., "An adaptive rational spectral method for differential
%   equations with rapidly varying solutions", Oxford DPhil Thesis, 2006.

if nargin < 2, k = 1; end

if N == 0, D = []; return, end
if N == 1, D = 0; return, end

% construct Chebyshev grid and weights
x = chebtech2.chebpts(N);
w = [.5 ; ones(N-1,1)]; w(2:2:end) = -1; w(N) = .5*w(N);

ii = (1:N+1:N^2)';              % indices of diagonal
Dx = bsxfun(@minus,x,x');       % all pairwise differences
Dx(ii) = Dx(ii) + 1;            % add identity
Dxi = 1./Dx;                    % reciprocal
Dw = bsxfun(@rdivide,w.',w);    % pairwise divisions
Dw(ii) = Dw(ii) - 1;            % subtract identity

% k = 1
D = Dw .* Dxi;
D(ii) = 0; D(ii) = - sum(D,2);              % negative sum trick

if k == 1, return, end

% k = 2
D = 2*D .* (repmat(D(ii),1,N) - Dxi);
D(ii) = 0; D(ii) = - sum(D,2);              % negative sum trick

% higher orders
for n = 3:k
    D = n*Dxi .* (Dw.*repmat(D(ii),1,N) - D);
    D(ii) = 0; D(ii) = - sum(D,2);          % negative sum trick
end

end
