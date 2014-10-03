function [x, w] = trigpts(n, dom)
%TRIGPTS    Equally spaced points.
%   TRIGPTS(N) returns N equispaced points in [-1, 1).
%
%   TRIGPTS(N, D), where D is vector of length 2 and N is a scalar integer,
%   returns N equispaced points in the interval [D(1),D(2)). 
%
%   [X, W] = TRIGPTS(N) or [X, W] = TRIGPTS(N, D) returns also a row vector 
%   of the weights for the trapezoidal rule.
%
% See also CHEBPTS, LEGPTS, JACPTS, LAGPTS, and HERMPTS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Special case (no points).
if ( nargin < 1 || n <= 0 )     
    x = []; 
    w = [];
    return
    
end

x = linspace(-pi, pi, n+1).';
x = x./pi;
x(end) = [];

if ( nargout > 1 )
    w = trigtech.quadwts(n);
end

if ( nargin > 1 )
    % Map to the right domain:
    x = diff(dom(1:2))*x/2 + mean(dom(1:2));
    
    if ( nargout > 1 )
        % Map to the right domain:
        w = w * diff(dom(1:2))/2;
    end
    
end

end