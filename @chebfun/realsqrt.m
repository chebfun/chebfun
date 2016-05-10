function R = realsqrt(X)
%REALSQRT   Real square root of a CHEBFUN.
%   REALSQRT(X) is the square root of the CHEBFUN of X.  An error is produced 
%   if X is negative or complex.
%
% See also CHEBFUN/SQRT, CHEBFUN/REALLOG, REALSQRT.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check for complex X:
if ( ~isreal(X) )
    error('CHEBFUN:CHEBFUN:realsqrt:complexRes', ...
        'REALSQRT produced complex result.');
end

% X is real positive, so call SQRT.
R = sqrt(X);

if ( normest(imag(R)) < 100*max(eps*vscale(R)) )
    R = real(R);
end

% Check for complex R:
if ( ~isreal(R) )
    error('CHEBFUN:CHEBFUN:realsqrt:complexRes', ...
        'REALSQRT produced complex result.');
end

end
