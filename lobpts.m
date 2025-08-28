function [x, w, v] = lobpts(n, alp, bet)
%LOBPTS   Gauss-Jacobi-Lobatto quadrature nodes and weights.
%  LOBPTS(N) returns N Legendre-Lobatto points X in [-1,1].
%
%  [X, W] = LOBPTS(N) returns also a row vector W of weights for
%  Gauss-Legendre-Lobatto quadrature.
%
%  [X, W, V] = LOBPTS(N) returns additionally a column vector V of weights in
%  the barycentric formula corresponding to the points X. The weights are
%  scaled so that max(abs(V)) = 1.
%
%  [...] = LOBPTS(N, ALP, BET) is similar, but for the Gauss-Jacobi-Lobatto
%  nodes and weights. Here ALP and BET should be scalars > -1.
%
%  In each case, N should be an integer greater than or equal to 2.
%
% See also CHEBPTS, LEGPTS, JACPTS, LEGPOLY, JACPOLY, RADAUPTS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer note:
%   The approach used here is to observe that the Gauss-Lobatto points are
%   precisely the roots of (1-x^2)P'_{n-1}(x), and that the roots of P'_{n-1}(x)
%   are the same as the roots of P^(1,1)_{n-2}(x) [NIST, (18.9.15)], which can
%   be obtained for JACPTS. A similar identity [NIST, (18.9.16)] is used for the
%   computation of the quadrature weights from those of JACPTS, and the missing
%   barycentric weights are determined by enforcing the interpolation of f(x) =
%   x and x^2 at x = 0.
%
%    x_j = roots of (1-x^2)P'_{n-1}(x)
%    w_j = { 2/(n*(n-1))                        : x_j = -1, 1
%          { 2/(n*(n-1)) * 1/[P_{n-1}(x_j)]^2   : otherwise
%
%   (Note that the weights for n-2 point Gauss-Jacobi with a = b = 1 satisfy 
%    u_j = C/(1-x_j^2)/[d/dx P^(1,1)_{n-2}(x_j)]^2)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO: Scaled domains?

if ( nargin < 3 )
    % Default to Gauss-Legendre-Lobatto:
    alp = 0; bet = 0;
end

%% Trivial cases:
if ( n == 1 )
    error('CHEBFUN:lobpts:notSupported', 'N = 1 is not supported.');
elseif ( n == 2 )
    x = [-1 ; 1];
    w = 2^(1+alp+bet)*[beta(bet+1, alp+2), beta(alp+1, bet+2)];
    v = [-1 ; 1];
    return
end

%% Call JACPTS():
[xi, w, v] = jacpts(n-2, alp+1, bet+1);

%% Nodes:
x = [-1 ; xi ; 1];

%% Quadrature weights:
w = [-1, w,  1];
w = w./(1-x.^2).';
if ( alp == 0 && bet == 0 )
    w([1 n]) = 2/(n*(n - 1));
else
    % The weights for the Gauss-Jacobi-Lobatto case are given explicitly by
    % Walter Gautschi, "High-order Gauss–Lobatto formulae", Numerical
    % Algorithms (2000).
    w(1) = 2^(1+alp+bet)*beta(bet+1, alp+n)*beta(bet+2,n-2)*(n-2);
    w(n) = 2^(1+alp+bet)*beta(alp+1, bet+n)*beta(alp+2,n-2)*(n-2);
end

%% Barycentric weights:
v = v./(1 - xi.^2);
v = v/max(abs(v));
rhs = -[sum(v) ; sum(v.*xi)];
vend = [1 1 ; -1 1]\rhs;
v = [vend(1) ; v ; vend(2)];

end