function I = integral(f, a, b, varargin)
%INTEGRAL   Evaluate integral of a CHEBFUN.
%   INTEGRAL(F, A, B) evaluates the integral of the CHEBFUN F over the interval
%   [A,B] using CHEBFUN/SUM(). If A and B are not given, the integral is
%   computed over the domain of F.
%
%   This function is a wrapper for CHEBFUN/SUM(). To use the original INTEGRAL()
%   on a CHEBFUN object, you can bypass this overloaded function by wrapping it
%   in an anonymous function:
%       integral(@(x) f(x), a, b);
%
% See also SUM.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin > 3 )
    % We'll allow this to slide...
    warning('CHEBFUN:CHEBFUN:integral:nargin', 'Too many input arguments.');
end

% Compute the integral with SUM:
if ( nargin == 1 )
    I = sum(f);
else
    I = sum(f, a, b);
end

end
