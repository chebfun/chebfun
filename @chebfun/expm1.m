function g = expm1(f)
%EXPM1   Compute EXP(X)-1 of a chebfun accurately.
%   EXPM1(Z) computes exp(Z)-1 accurately in the case where the chebfun Z is
%   small on its domain. Complex Z is accepted.
%
% See also EXPM1, CHEBFUN/LOG1P.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfun.pref();
end

% Call the compose method:
g = compose(f, @expm1, pref);

end