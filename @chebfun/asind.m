function F = asind(F, pref)
%ASIND   Inverse sine of a CHEBFUN, result in degrees.
%   ASIND(F) computes the inverse sine (in degrees) of the CHEBFUN F.
%
%   ASIND(F, PREF) does the same but uses the CHEBPREF object PREF when
%   computing the composition.
%
% See also SIND, ASIN.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebpref();
end

% Loop over the columns of F:
for k = 1:numel(F)
    % Call the compose method:
    F(k) = compose(F(k), @asind, pref);
end

end
