function g = tanh(f, pref)
%TANH   Hyperbolic tangent of a CHEBFUN.
%   TANH(F) computes the hyperbolic tangent of the CHEBFUN F.
%
%   TANH(F, PREF) does the same but uses the CHEBPREF object PREF when
%   computing the composition.
%
% See also ATANH.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebpref();
end

% Call the compose method:
g = compose(f, @tanh, pref);

end