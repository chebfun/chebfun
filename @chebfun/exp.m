function F = exp(F, varargin)
%EXP   Exponential of a CHEBFUN.
%   EXP(F) computes the exponential of the CHEBFUN F.
%
%   EXP(F, PREF) does the same but uses the CHEBFUNPREF object PREF when computing
%   the composition.
%
% See also EXPM1.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call the compose method:
F = compose(F, @exp, varargin{:});

end
