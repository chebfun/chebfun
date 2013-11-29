function [x, w, v] = radaupts(n, varargin)
%RADAUPTS   Gauss-Legendre-Radau Quadrature Nodes and Weights.
%  RADAUPTS(N) returns N Legendre-Radau points X in (-1,1).
%
%  [X, W] = RADAUPTS(N) returns also a row vector W of weights for
%  Gauss-Legendre-Lobatto quadrature.
%
%  [X, W, V] = RADAUPTS(N) returns additionally a column vector V of weights in
%  the barycentric formula corresponding to the points X. The weights are scaled
%  so that max(abs(V)) = 1.
%
%  In each case, N should be a positive integer.
%
%  See also CHEBPTS, LEGPTS, JACPTS, LEGPOLY, LOBPTS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer note:
%   The approach used here is to observe that the Gauss-Radau points are
%   precisely the roots of (1+x)P^(0,1)_{n-2}(x), which can be obtained from
%   JACPTS. A similar identity is used for the computation of the quadrature
%   weights from those of JACPTS, and the missing barycentric weights are
%   determined by enforcing the interpolation of x at x = 0.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Trivial cases:
if ( n == 1 )
    x = -1;
    w = 2;
    v = 1;
    return
elseif ( n == 1 )
    x = [-1, 1/3];
    w = [.5 ; 1.5];
    v = [-1 ; 1];
    return
end

%% Call JACPTS():
[x, w, v] = jacpts(n - 1, 0, 1, varargin{:});

%% Nodes:
x = [-1 ; x];

%% Quadrature weights:
w = [2/n^2, w./(1 + x(2:end).')];

%% Barycentric weights:
v = v./(1 + x(2:end));
v = v/max(abs(v));
v1 = -abs(sum(x(2:end).*v));
v = [v1 ; v];

end

