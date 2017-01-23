function f = lagrange(x, varargin)
%LAGRANGE   Compute Lagrange basis functions.
%   F = CHEBFUN.LAGRANGE(X) returns a CHEBFUN object F representing the Lagrange
%   polynomials for the points X(0), ..., X(N). That is, each column of F is a
%   a polynomial of degree N which satisfies F(X,:) = eye(length(X)).
%
%   F = CHEBFUN.LAGRANGE(X, DOM) restricts the result F to the domain DOM. DOM
%   _must_ be passed if X is a scalar.
%
% See also INTERP1, VANDER.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check inputs:
n = length(x);
if ( n == 0 )
    f = chebfun();
    return
elseif ( n == 1 && nargin < 2 )
    error('CHEBFUN:CHEBFUN:lagrange:nodomain', ...
        'Domain must be specfied when X is a scalar.')
end

% Check for uniqueness:
if ( length(unique(x)) ~= n )
    error('CHEBFUN:CHEBFUN:lagrange:nonunique', ...
        'Interpolation points must be unique.')
end

% Make interpolation data (identity matrix):
y = eye(n); 

% X values must be sorted for INTERP1:
[x, idx] = sort(x);
y = y(:,idx);

% Call INTERP1():
f = chebfun.interp1(x, y, 'poly', varargin{:});

end
