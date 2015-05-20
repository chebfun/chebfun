function [x, w] = trigpts(n)
%TRIGPTS   Equispaced points in [-1, 1).
%   TRIGPTS(N) returns N equispaced points in [-1, 1).
%
%   [X, W] = TRIGPTS(N) returns also a row vector of the weights for
%   the trapezoidal rule.
%
% See also CHEBPTS, LEGPTS, JACPTS, LAGPTS, and HERMPTS.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Special case (no points).
if ( n <= 0 )     
    x = []; 
    w = [];  
    return
    
end

x = linspace(-1, 1, n+1).';
x(end) = [];

% Quadrature weights:
if ( nargout > 1 )
    w = trigtech.quadwts(n);
end
    
end
