function F = sec(F, pref)
%SEC   Secant of a CHEBFUN.
%   SEC(F) computes the secant of the CHEBFUN F.
%
%   SEC(F, PREF) does the same but uses the CHEBPREF object PREF when computing
%   the composition.
%
% See also ASEC, SECD.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebpref();
end

% Loop over the columns of F:
for k = 1:numel(F)
    % Call the compose method:
    F(k) = compose(F(k), @sec, pref);
end

end