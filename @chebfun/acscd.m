function F = acscd(F, pref)
%ACSCD   Inverse cosecant of a CHEBFUN, result in degrees.
%   ACSCD(F) computes the inverse cosecant (in degrees) of the CHEBFUN F.
%
%   ACSCD(F, PREF) does the same but uses the CHEBPREF object PREF when
%   computing the composition.
%
% See also CSCD, ACSC.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebpref();
end

% Loop over the columns of F:
for k = 1:numel(F)
    % Call the compose method:
    F(k) = compose(F(k), @acscd, pref);
end

end
