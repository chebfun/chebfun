function D = baryDiffMat(x, w, k, t)
%BARYDIFFMAT  Barycentric differentiation matrix.
%   D = BARYDIFFMAT(N, X) is the matrix that maps function values at the points
%   X to values of the derivative of the interpolating polynomial at those
%   points.
%
%   D = BARYDIFFMAT(X, W) is the same, but here the derivative is of the
%   barycentric interpolant defined by the points X and the weights W.
%
%   D = BARYDIFFMAT(X, W, K) is the same, but for the Kth derivative.
%
%   D = BARYDIFFMAT(X, W, K, T) but accepts T = ACOS(X), which allows more
%   accurate computation of pairwise differences of X in some cases. (See [4]).
%
%   The matrices are computed using the 'hybrid' formula of Schneider & Werner
%   [1] and Welfert [2] proposed by Tee [3] with the tricks suggested by [4].
%
% References:
%  [1] Schneider, C. and Werner, W., "Some new aspects of rational
%   interpolation", Math. Comp. (47) 285--299, 1986.
%  [2] Welfert, B. D., "Generation of pseudospectral matrices I", SINUM, (34)
%   1640--1657.
%  [3] Tee, T. W., "An adaptive rational spectral method for differential
%   equations with rapidly varying solutions", Oxford DPhil Thesis, 2006.
%  [4] R. Baltensperger and M.R. Trummer, "Spectral Differencing with a Twist",
%   SIAM J. Sci. Comp., Vol. 24, No. 5, pp. 1465–1487, 2003.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

N = length(x);

%% Parse inputs and check for trivial cases:
if ( N == 0 )
    D = [];
    return
elseif ( N == 1 )
    D = 0;
    return
end
if ( nargin < 3 )
    k = 1;
end
if ( k == 0 )
    D = eye(N);
    return
end
if ( nargin < 2 )
    w = baryWeights(x);
end

if ( length(x) ~= length(w) )
    error('CHEBFUN:COLLOC:baryDiffMat:xwlengths', ...
        'length(x) must equal length(w).');
end

%% Construct Dx and Dw:
ii = (1:N+1:N^2).';                                 % Indices of diagonal.
if ( nargin == 4 ) 
    % Trig identity as described in [4]:
    t = flipud(t);
    Dx = 2*bsxfun( @(t, tp) sin(t + tp).*sin(t - tp), t/2, t.'/2 ); 
else
    Dx = bsxfun(@minus, x, x.');                    % All pairwise differences.
end

% Flipping trick in [4]:
DxRot = rot90(Dx, 2);
idxTo = rot90(~triu(ones(N)));
Dx(idxTo) = -DxRot(idxTo);

Dx(ii) = 1;                                         % Add identity.
Dxi = 1./Dx;                                        % Reciprocal.
Dw = bsxfun(@rdivide, w.', w);                      % Pairwise divisions.
Dw(ii) = 0;                                         % Subtract identity.

%% k = 1
D = Dw .* Dxi;
D(ii) = 0; D(ii) = -sum(D, 2);                      % Negative sum trick.

% Forcing symmetry for even N:
D(ii(end:-1:N-floor(N/2)+1)) = -D(ii(1:floor(N/2)));

if ( k == 1 )
    return
end

%% k = 2
D = 2*D .* ( repmat(D(ii), 1, N) - Dxi );
D(ii) = 0; D(ii) = -sum(D, 2);                      % Negative sum trick.

%% k = 3...
for n = 3:k
    D = n*Dxi .* ( Dw.*repmat(D(ii), 1, N) - D );
    D(ii) = 0; D(ii) = -sum(D, 2);                  % Negative sum trick.
end

end
