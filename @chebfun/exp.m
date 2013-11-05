function g = exp(f, pref)
%EXP   Exponential of a CHEBFUN.
%   EXP(F) computes the exponential of the CHEBFUN F.
%
%   EXP(F, PREF) does the same but uses the CHEBPREF object PREF when
%   computing the composition.
%
%   See also EXPM1.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebpref();
end

% Call the compose method:
g = compose(f, @exp, pref);

end
