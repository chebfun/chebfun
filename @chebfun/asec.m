function F = asec(F, varargin)
%ASEC   Inverse secant of a CHEBFUN.
%   ASEC(F) computes the inverse secant of the CHEBFUN F.
%
%   ASEC(F, PREF) does the same but uses the CHEBFUNPREF object PREF when
%   computing the composition.
%
% See also SEC, ASECD.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% Call the compose method:
F = compose(F, @asec, varargin{:});

end
