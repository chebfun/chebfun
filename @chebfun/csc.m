function F = csc(F, pref)
%CSC   Cosecant of a CHEBFUN.
%   CSC(F) computes the cosecant of the CHEBFUN F.
%
%   CSC(F, PREF) does the same but uses the CHEBPREF object PREF when
%   computing the composition.
%
% See also ACSC, CSCD.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebpref();
end

% Loop over the columns of F:
for k = 1:numel(F)
    % Call the compose method:
    F(k) = compose(F(k), @cos, pref);
end

end