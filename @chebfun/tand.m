function F = tand(F, pref)
%TAND   Tangent of a CHEBFUN, result in degrees.
%   TAND(F) computes the tangent (in degrees) of the CHEBFUN F.
%
%   TAND(F, PREF) does the same but uses the CHEBPREF object PREF when
%   computing the composition.
%
% See also ATAND, TAN.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebpref();
end

% Loop over the columns of F:
for k = 1:numel(F)
    % Call the compose method:
    F(k) = compose(F(k), @tand, pref);
end

end
