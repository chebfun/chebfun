function g = sinc(f, pref)
%SINC   Sinc function of a chebfun.
%
% See also SIN.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfun.pref();
end

% Call the compose method:
g = compose(f, @(x) sin(pi*x)./(pi*x), pref);

end
