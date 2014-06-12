function x = fourpts(n, dom)
%FOURPTS    Fourier points.
%   FOURPTS(N) returns N equispaced points in [-1, 1).
%
%   FOURPTS(N, D), where D is vector of length 2 and N is a scalar integer,
%   returns N equispaced points in the interval [D(1),D(2)). 
%
% See also CHEBPTS, LEGPTS, JACPTS, LAGPTS, and HERMPTS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

x = linspace(-pi, pi, n+1).';
x = x./pi;
x(end) = [];

% Parse inputs:
if ( nargin == 1 )
    dom = [-1, 1]; 
end

% Map to the right domain: 
x = diff(dom(1:2))*x/2 + mean(dom(1:2)); 
    
end