function g = cscd(f, pref)
%CSCD   Cosecant of a CHEBFUN, result in degrees.
%   CSCD(F) computes the cosecant (in degrees) of the CHEBFUN F.
%
%   CSCD(F, PREF) does the same but uses the preference structure PREF when
%   computing the composition.
%
% See also ACSCD, CSC.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfun.pref();
end

% Call the compose method:
g = compose(f, @cscd, pref);

end
