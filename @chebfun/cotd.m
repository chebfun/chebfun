function g = cotd(f, pref)
%COSD   Cotangent of a chebfun, result in degrees.
%
% See also ACOTD, COT.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfun.pref();
end

% if ( ~isfinite(f) )
%     error('CHEBFUN:cotd:inf',...
%         'COTD is not defined for functions which diverge to infinity');
% end

% Call the compose method:
g = compose(f, @cotd, pref);

end