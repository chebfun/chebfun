function [x, w, v] = lobpts(n, varargin)
%LOBPTS   Gauss-Legendre-Lobatto Quadrature Nodes and Weights.
%  LOBPTS(N) returns N Legendre-Lobatto points X in [-1,1].
%
%  [X, W] = LOBPTS(N) returns also a row vector W of weights for
%  Gauss-Legendre-Lobatto quadrature.
%
%  [X, W, V] = LOBPTS(N) returns additionally a column vector V of weights in
%  the barycentric formula corresponding to the points X. The weights are
%  scaled so that max(abs(V)) = 1.
%
%  In each case, N should be an integer greater than or equal to 2.
%
% See also CHEBPTS, LEGPTS, JACPTS, LEGPOLY, RADAUPTS.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer note:
%   The approach used here is to observe that the Gauss-Lobatto points are
%   precisely the roots of (1-x^2)P'_{n-1}(x), and that the roots of P'_{n-1}(x)
%   are the same as the roots of P^(1,1)_{n-2}(x) [NIST, (18.9.15)], which can
%   be obtained for JACPTS. A similar identity [NIST, (18.9.16)] is used for the
%   computation of the quadrature weights from those of JACPTS, and the missing
%   barycentric weights are determined by enforcing the interpolation of f(x) =
%   x or x^2 at x = 0 in the even or odd case respectively.
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

%% Trivial cases:
if ( n == 1 )
    error('CHEBFUN:lobpts:notSupported', 'N = 1 is not supported.');
elseif ( n == 2 )
    x = [-1 ; 1];
    w = [1, 1];
    v = [-1 ; 1];
    return
elseif ( n == 3 )
    x = [-1 ; 0 ; 1];
    w = [1, 4, 1]/3;
    v = [-.5 ; 1 ; -.5];
    return
end

%% Call JACPTS():
[x, w, v] = jacpts(n - 2, 1, 1, varargin{:});

%% Nodes:
x = [-1 ; x ; 1];

%% Quadrature weights:
w = [-1, w,  1];
w = w./(1-x.^2).';
w([1 end]) = 2/(n*(n - 1));

%% Barycentric weights:
v = v./(1 - x(2:n-1).^2);
v = v/max(abs(v));
if ( mod(n, 2) )
    v1 = -abs(sum(v.*x(2:end-1).^2)/2);
    sgn = 1;
else
    v1 = -abs(sum(v.*x(2:end-1))/2);
    sgn = -1;
end
v = [v1 ; v ; sgn*v1];

end

