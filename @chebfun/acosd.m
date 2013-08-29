function g = acosd(f, pref)
%ACOSD   Cosine of a CHEBFUN, result in degrees.
%   ACOSD(F) computes the cosine (in degrees) of the CHEBFUN F.
%
%   ACOSD(F, PREF) does the same but uses the preference structure PREF when
%   computing the composition.
%
% See also ACOS, COS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfun.pref();
end

% Call the compose method:
g = compose(f, @acosd, pref);

end
