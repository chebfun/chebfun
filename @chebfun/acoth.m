function F = acoth(F, pref)
%ACOTH   Inverse hyperbolic cotangent of a CHEBFUN.
%   ACOTH(F) computes the inverse hyperbolic cotangent of the CHEBFUN F.
%
%   ACOTH(F, PREF) does the same but uses the CHEBPREF object PREF when
%   computing the composition.
%
% See also COTH.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebpref();
end

% Loop over the columns of F:
for k = 1:numel(F)
    % Call the compose method:
    F(k) = compose(F(k), @acoth, pref);
end

end
