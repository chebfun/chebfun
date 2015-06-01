function w = wronskian(L, varargin)
%WRONSKIAN   Wronskian of chebfuns.
%   WRONSKIAN(L, f1, ..., fn) computes the wronskian of the CHEBFUN objects fk.
%   L is the linear differential operator and f1, ..., fn are solutions of the
%   homogenous problem. The algorithm is based on Abel's identity for computing
%   the wronskian.
%   
%   WRONSKIAN(L, F) does the same where F is a quasimatrix or an array-valued
%   CHEBFUN or a CHEBMATRIX. If F is a CHEBMATRIX, then the form of F is assumed
%   to be a CHEBMATRIX with a single column with n CHEBFUNs.
%
% See also LINOP

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% At the moment, we do not support systems of ODEs:
if ( min(size(L.blocks)) > 1 )
    error('CHEBFUN:LINOP:wronskian:system', 'ODE systems not supported.');
end

if ( nargin > 2 )
    % Input is a comma separated list of CHEBFUNs.
    F = horzcat(varargin{:});
else
    % F is either a quasimatrx, array-valued Chebfun, or chebmatrix.
    F = varargin{1};
end

% Extract the coefficients of the operator:
c = toCoeff(L.blocks{1});
% The (n-1)th coefficient in standard form:
p = c{2}./c{1};
dom = L.domain;
a = dom(1);
n = size(c, 2) - 1;
nFuns = size(F, 2);

if ( nFuns ~= n )
    error('CHEBFUN:LINOP:wronskian:nFuns', ...
        'Number of functions is not the same as the order of the operator.')
end

W = zeros(n);
for i = 1:n
    W(i,:) = feval(F, a);
    F = diff(F); 
end

% Compute the determinant at the left end of the domain:
A = det(W);

% Apply Abel's identity:
w = A*exp(-cumsum(p));

% Ensure w(a) = A:
w = (w - w(a)) + A;

end
    
    
