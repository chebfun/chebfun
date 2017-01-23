function f = volt(K, v)
%VOLT      Volterra integral operator for CHEBFUN2.
%   V = VOLT(K, f) returns a row CHEBFUN resulting from the integral
%
%      f(x) = (K*v)(x) = int( K(x,y) v(y), y=a..x ),
%
%   where K is defined on a domain [a,b]x[a,b].
%
%   The kernel function K(x,y) must be a smooth CHEBFUN2 defined on a square
%   domain.
%
% Example:
%   f = volt(chebfun2(@(x,y) exp(x-y)),chebfun('x'));
%
% See also FRED.

% Copyright 2017 by The University of Oxford and The Chebfun2 Developers.
% See http://www.chebfun.org/ for Chebfun2 information.

if ( ~isa(K, 'chebfun2') )
    error('CHEBFUN:CHEBFUN2:volt:input', ...
          'First argument must be a CHEBFUN2.');
end

% Get the low rank representation for f.
[cols, D, rows] = cdr(K);
dom = K.domain;

if isa(v, 'function_handle')
    % Convert to a CHEBFUN on the right interval.
    v = chebfun(v, dom(3:4));
end

% Domain compatibility:
if ( ~domainCheck(cols, v) )
    error('CHEBFUN:CHEBFUN2:volt:domainMismatch', ...
          'Domain of CHEBFUN and CHEBFUN2 kernel do not match.');
end

RR = rows * D;

% Cumsum with cols and v:  (This can be sped up.)
f = chebfun(0, dom(1:2));
for jj = length(K) : -1 : 1
    CC = cumsum(v .* cols(:,jj));
    f = f + CC .* RR(:,jj);
end

% Transpose to preserve linear algebra:
f = f.';

end
