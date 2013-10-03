function g = sinc(f, pref)
%SINC   Sinc function of a CHEBFUN.
%   SINC(F) computes the sinc function of the CHEBFUN F.
%
%   SINC(F, PREF) does the same but uses the preference structure PREF when
%   computing the composition.
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
