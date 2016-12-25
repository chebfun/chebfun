function f = trigpoly(n, d)
%TRIGPOLY   Trigonometric polynomial in complex exponential form.
%   F = TRIGPOLY(N) returns the degree N trigonometric polynomial exp(1i*pi*N*x) 
%   on [-1,1], where N may be a vector of integers.
%
%   F = TRIGPOLY(N, D), where D is an interval or a domain, gives the same
%   result scaled accordingly.
%
% See also CHEBPOLY, LEGPOLY, and TRIGPTS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information. 

% Parse input
if ( any(mod(n, 1) ~= 0) )
    error('CHEBFUN:trigpoly:integern', ...
    'The first argument must be a vector of integers.');
end

if ( nargin == 1 )
    d = chebfunpref().domain;
end    

% Cannot handle unbounded domains:
if ( any(isinf(d)) )
    error('CHEBFUN:trigpoly:infdomain', ...
    'Trigonometric polynomial construction is not allowed on unbounded domain.');
end

% Cannot handle interior breakpoints:
if ( numel(d) > 2 )
    error('CHEBFUN:trigpoly:breakpoints', ...
    'Trigonometric polynomials can not be constructed on domains with break points.');
end

% Construct the trigonometric Fourier coefficients:
N = max(abs(n))+1;
c = zeros(2*N-1, numel(n));
for k = 1:numel(n)
    c(N-n(k), k) = 1;
end

% Construct a CHEBFUN from the coefficients:
f = chebfun(c, d([1, end]), 'periodic', 'coeffs');

% Transpose if required:
if ( size(n, 1) > 1 )
    f = f.'; 
end

end

