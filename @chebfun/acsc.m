function F = acsc(F, pref)
%ACSC   Inverse cosecant of a CHEBFUN.
%   ACSC(F) computes the inverse cosecant of the CHEBFUN F.
%
%   ACSC(F, PREF) does the same but uses the CHEBPREF object PREF when
%   computing the composition.
%
% See also CSC, ACSCD.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebpref();
end

% Call the compose method:
F = compose(F, @acsc, pref);

end
