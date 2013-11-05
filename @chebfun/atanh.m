function g = atanh(f, pref)
%ATANH   Inverse hyperbolic tangent of a CHEBFUN.
%   ATANH(F) computes the inverse hyperbolic tangent of the CHEBFUN F.
%
%   ATANH(F, PREF) does the same but uses the CHEBPREF object PREF when
%   computing the composition.
%
% See also TANH.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebpref();
end

% Call the compose method:
g = compose(f, @atanh, pref);

end
