function [x, w, v] = radaupts(n, alp, bet, dom)
%RADAUPTS   Gauss-Jacobi-Radau quadrature nodes and weights.
%  RADAUPTS(N) returns N Legendre-Radau points X in [-1,1).
%
%  [X, W] = RADAUPTS(N) returns also a row vector W of weights for
%  Gauss-Legendre-Radau quadrature.
%
%  [X, W, V] = RADAUPTS(N) returns additionally a column vector V of weights in
%  the barycentric formula corresponding to the points X. The weights are scaled
%  so that max(abs(V)) = 1.
%
%  [...] = RADUAPTS(N, ALP, BET) is similar, but for the Gauss-Jacobi-Radau
%  nodes and weights. Here ALP and BET should be scalars > -1.
% 
%  [...] = RADAUPTS(N, [A,B]) or [...] = RADAUPTS(N, ALP, BET, [A,B])
%  scales the nodes and weights for the interval [A,B) .
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

% Parse inputs:
if ( nargin < 2 )
    % Default to Gauss-Legendre-Radau on [-1 1]:
    alp = 0; bet = 0; dom = [-1 1];
elseif ( nargin == 2 && numel(alp) == 2 )
    % Default to Gauss-Legendre-Radau on [A, B]:
    dom = alp; alp = 0; bet = 0;
elseif ( nargin == 3 )
    % Default to Gauss-Jacobi-Radau on [-1,1]:
    dom = [-1 1];
elseif ( nargin ~= 4 )
    error('CHEBFUN:radaupts:inputs', 'Unable to parse inputs.')
end

% Check inputs:
if ( ~isscalar(alp) ||  alp <= -1 )
    error('CHEBFUN:radaupts:alp', 'ALP > 1 must be a scalar.');
end
if ( ~isscalar(bet) ||  bet <= -1 )
    error('CHEBFUN:radaupts:bet', 'BET > 1 must be a scalar.');
end
if ( numel(dom)~=2 )
    error('CHEBFUN:radaupts:domain', 'Invalid domain argument.');
end

%% Trivial cases:
if ( n == 1 )
    x = -1;
    w = 2^(1+alp+bet)*beta(1+alp,1+bet);
    v = 1;
else

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

% Scale the nodes and weights:
if ( dom(1) == -1 && dom(2) == 1 )
    % Nodes are already on [-1, 1];
else
    x = dom(2)*(x + 1)/2 + dom(1)*(1 - x)/2;
    w = (diff(dom)/2)*w;
end

end

