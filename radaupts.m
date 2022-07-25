function [x, w, v] = radaupts(n, alp, bet)
%RADAUPTS   Gauss-Jacobi-Radau quadrature nodes and weights.
%  RADAUPTS(N) returns N Legendre-Radau points X in [-1,1).
%
%  [X, W] = RADAUPTS(N) returns also a row vector W of weights for
%  Gauss-Legendre-Lobatto quadrature.
%
%  [X, W, V] = RADAUPTS(N) returns additionally a column vector V of weights in
%  the barycentric formula corresponding to the points X. The weights are scaled
%  so that max(abs(V)) = 1.
%
%  [...] = RADUAPTS(N, ALP, BET) is similar, but for the Gauss-Jacobi-Radau
%  nodes and weights. Here ALP and BET should be scalars > -1.
%
%  In each case, N should be a positive integer.
%
% See also CHEBPTS, LEGPTS, JACPTS, LEGPOLY, JACPOLY, LOBPTS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer note:
%   The approach used here is to observe that the Gauss-Radau points are
%   precisely the roots of (1+x)P^(0,1)_{n-2}(x), which can be obtained from
%   JACPTS. A similar identity [NIST, (18.9.5) and (18.9.17)] is used for the
%   computation of the quadrature weights from those of JACPTS, and the missing
%   barycentric weights are determined by enforcing the interpolation of f(x) =
%   x at x = 0.
%
%    x_j = roots of (1+x)P^(0,1)_{n-2}(x)
%    w_j = { 2/n^2                                 : x_j = -1
%          { 1/(1-x_j) * 1/[d/dx P_{n-1}(x_j)]^2   : otherwise
%
%   (Note that the weights for n-1 point Gauss-Jacobi with a = 0, b = 1 satisfy 
%    u_j = C/(1-x_j^2)/[d/dx P^(0,1)_{n-1}(x_j)]^2)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO: Scaled domains?

if ( nargin < 3 )
    % Default to Gauss-Legendre-Lobatto:
    alp = 0; bet = 0;
end

%% Trivial cases:
if ( n == 1 )
    x = -1;
    w = 2^(1+alp+bet)*beta(1+alp,1+bet);
    v = 1;
    return
end

%% Call JACPTS():
[xi, w, v] = jacpts(n-1, alp, bet+1);

%% Nodes:
x = [-1 ; xi];

%% Quadrature weights:
wi = w./(1 + xi.');
if ( alp == 0 && bet == 0 )
    w = [2/n^2, wi];
else
    % See Walter Gautschi, "Gauss–Radau formulae for Jacobi and Laguerre
    % weight functions", Mathematics and Computers in Simulation, (2000).
    w = [2^(alp+bet+1)*beta(bet+1,n)*beta(alp+n,bet+1)*(bet+1), wi];
end

%% Barycentric weights:
v = v./(1 + xi);
v = v/max(abs(v));
v1 = -sum(v);
v = [v1 ; v];

end

