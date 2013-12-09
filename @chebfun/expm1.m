function F = expm1(F, pref)
%EXPM1   Compute EXP(F)-1 of a CHEBFUN accurately.
%   EXPM1(F) computes EXP(F)-1 accurately in the case where the CHEBFUN F is
%   small on its domain. Complex F is accepted.
%
%   EXPM1(F, PREF) does the same but uses the CHEBPREF object PREF when
%   computing the composition.
%
% See also EXP, LOG1P.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebpref();
end

% Loop over the columns of F:
for k = 1:numel(F)
    % Call the compose method:
    F(k) = compose(F(k), @expm1, pref);
end

end
