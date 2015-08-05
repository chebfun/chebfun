function F = asecd(F, varargin)
%ASECD   Inverse secant of a CHEBFUN, result in degrees.
%   ASECD(F) computes the inverse secant (in degrees) of the CHEBFUN F.
%
%   ASECD(F, PREF) does the same but uses the CHEBFUNPREF object PREF when
%   computing the composition.
%
% See also SECD, ASEC.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% Call the compose method:
F = compose(F, @asecd, varargin{:});

end
