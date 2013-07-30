function I = integral(f, a, b, varargin)
%INTEGRAL  Numerically evaluate integral of a chebfun.
%   INTEGRAL(F, A, B) evaluates the integral of the chebfun F over the interval
%   [A,B] using chebfun/sum. If A and B are not given, the integral is computed
%   over the domain of F.
%
%   This function is a wrapper for CHEBFUN/SUM(). To use the original INTEGRAL()
%   on a chebfun, you can bypass this overloaded function by wrapping it in an
%   anonymous function:
%       integral(@(x) f(x), a , b);
%
%   See also INTEGRAL, CHEBFUN/SUM

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin > 3 )
    % We'll allow this to slide...
    warning('CHEBFUN:integral:nargin', 'Too many input arguments.');
end

% Compute the integral with SUM:
if ( nargin == 1 )
    I = sum(f);
else
    I = sum(f, a, b);
end
