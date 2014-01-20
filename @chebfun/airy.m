function F = airy(K, F, scale)
%AIRY   Airy function of a CHEBFUN.
%   AIRY(F) returns the Airy function Ai(F) of a CHEBFUN F.
%
%   AIRY(K, F) returns various Airy functions specified by K:
%     0 - (default) is the same as airy(Z)
%     1 - returns the derivative, Ai'(Z)
%     2 - returns the Airy function of the second kind, Bi(Z)
%     3 - returns the derivative, Bi'(Z)
%
%   AIRY(K, F, SCALE) returns a scaled AIRY(K, F) specified by SCALE:
%     0 - (default) is that same as AIRY(K, Z)
%     1 - returns airy(K, F) scaled by EXP(2/3*F.^(3/2)) for K = 0, 1,
%         and scaled by EXP(-ABS(2/3.*REAL(F.^(3/2)))) for K = 2, 3.
%
% See also BESSELJ.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% TODO: Preference input.

% Parse the inputs as described in help text:
if ( nargin == 1 )
    F = K;
    K = 0;
end
if ( nargin < 3 )
    scale = 0;
end

for k = 1:numel(F)
    F(k) = columnAiry(K, F(k), scale);
end

end

function g = columnAiry(K, f, scale)

% The standard Airy function:
g = compose(f, @(x) airy(K, x));

if ( scale == 0 )
    % Standard case (no scaling).

elseif ( (scale == 1) && ((K == 0) || (K == 1)) )
    % Scaled with k = 1, 2:
    scl = exp(2/3*f.^(3/2));
    g = scl.*g;

elseif ( (scale == 1) && ((K == 2) || (K == 3)) )
    % Scaled with k = 2, 3:
    scl = exp(-abs(2/3*real(f.^(3/2))));
    g = scl.*g;

else
    % Invalid parameter sequence:
    error('CHEBFUN:airy:params', 'Invalid paramter selection.');

end

end
