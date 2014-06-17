function F = expm1(F, varargin)
%EXPM1   Compute EXP(F)-1 of a CHEBFUN accurately.
%   EXPM1(F) computes EXP(F)-1 accurately in the case where the CHEBFUN F is
%   small on its domain. Complex F is accepted.
%
%   EXPM1(F, PREF) does the same but uses the CHEBFUNPREF object PREF when
%   computing the composition.
%
% See also EXP, LOG1P.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call the compose method:
F = compose(F, @expm1, varargin{:});

end
