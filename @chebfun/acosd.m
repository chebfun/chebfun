function F = acosd(F, pref)
%ACOSD   Cosine of a CHEBFUN, result in degrees.
%   ACOSD(F) computes the cosine (in degrees) of the CHEBFUN F.
%
%   ACOSD(F, PREF) does the same but uses the CHEBPREF object PREF when
%   computing the composition.
%
% See also ACOS, COS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebpref();
end

% Loop over the columns of F:
for k = 1:numel(F)
    % Call the compose method:
    F(k) = compose(F(k), @acosd, pref);
end

end
