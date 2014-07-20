function w = wronskian(L, varargin)
%WRONSKIAN   Wronskian of chebfuns.
%   WRONSKIAN(L, f1, ..., fn) computes the wronskian of the chebfuns.
%   L is the linear differential operator and f1, ..., fn are solutions of the
%   homogenous problem. The algorithm is based on Abel's identity for computing
%   the wronskian.
%
% See also LINOP

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Extract the coefficients of the operator:
L = linop(L);
c = toCoeff(L.blocks{1});
% The (n-1)th coefficient in standard form:
p = c{2}./c{1};
dom = L.domain;
a = dom(1);
sz = size(c);
n = sz(2)-1;
if ( length(varargin) ~= n )
    error( 'CHEBFUN:WRONSKIAN', 'Number of chebfuns is not the same as the order of the operator' )
end

% [TODO]: This can be vectorized using fancy chembatrices etc?
W = zeros(n);
for i = 1:n
    for j = 1:n
        f = varargin{j};
        W(i, j) = f(a);
        varargin{j} = diff(f);
    end
end

% Compute the determinant at the left end of the domain:
A = det(W);

% Apply Able's identity
w = A*exp(-cumsum(p));
% Make sure w(a) = A
w = (w - w(a)) + A;
end
    
    