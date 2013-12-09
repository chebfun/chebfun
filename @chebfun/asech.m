function F = asech(F, pref)
%ASECH   Inverse hyperbolic secant of a CHEBFUN.
%   ASECH(F) computes the inverse hyperbolic secant of the CHEBFUN F.
%
%   ASECH(F, PREF) does the same but uses the CHEBPREF object PREF when
%   computing the composition.
%
% See also SECH.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebpref();
end

% Loop over the columns of F:
for k = 1:numel(F)
    % Call the compose method:
    F(k) = compose(F(k), @asech, pref);
end

end
