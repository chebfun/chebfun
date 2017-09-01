function K = intmat(N, varargin)
%INTMAT   Spectral integration matrix.
%   K = INTMAT(N) returns the N x N integration matrix associated with the
%   Chebyshev spectral collocation method at second-kind Chebyshev points. 
%
%   K = INTMAT(N, P) returns the N x N integration matrix of order P.
%
%   K = INTMAT(N, P, DOM) scales the integration matrix K to the domain
%   DOM. DOM should be a 1x2 vector.
%
%   K = INTMAT([M N]) returns the M x N first-order rectangular integration 
%   matrix which maps from an N-point Chebyshev grid of the second kind to an 
%   M-point Chebyshev grid of the same kind.
%   
%   K = INTMAT([M N], P) returns an M x N rectangular integration matrix of 
%   order P which maps from an N-point to an M-point Chebyshev grid, both of
%   second kind.
%
%   K = INTMAT([M N], P, DOM) returns the same K but scaled to the domain DOM.
%
% See also CUMSUMMAT, DIFFMAT, DIFFROW, INTROW.

% This code has been written for Aurentz and Trefethen, 
% "Block operators and spectral discretizations", and is not
% yet fully compliant with standard Chebfun coding practices.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse inputs
[m, n, p, dom] = parseInputs(N, varargin{:});

% Build Lagrange basis
K = chebfun(eye(n),dom); 

% Integrate p times
K = cumsum(K,p);

% Evaluate at grid of size m
K = feval(K,chebpts(m,dom));

end

function [m, n, p, dom] = parseInputs(N, varargin)
% Parse the inputs to INTMAT.

p = 1;
dom = [-1 1];

if ( isscalar(N) )
    n = N;
    m = n;
else
    m = N(1);
    n = N(2);
end

for j = 1:numel(varargin)
    v = varargin{j};
    if ( isnumeric(v) )
        if ( isscalar(v) )
            p = v;
        else
            dom = v;
        end
    else
        error('CHEBFUN:intmat:unknown', ...
            'Unknown input of type %s.', class(v));
    end
end

end

