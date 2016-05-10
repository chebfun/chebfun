function F = sec(F, varargin)
%SEC   Secant of a CHEBFUN.
%   SEC(F) computes the secant of the CHEBFUN F.
%
%   SEC(F, PREF) does the same but uses the CHEBFUNPREF object PREF when computing
%   the composition.
%
% See also ASEC, SECD.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call the compose method:
F = compose(F, @sec, varargin{:});

end
