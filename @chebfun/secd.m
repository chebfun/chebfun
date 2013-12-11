function F = secd(F, varargin)
%SECD   Secant of a CHEBFUN, result in degrees.
%   SECD(F) computes the secant (in degrees) of the CHEBFUN F.
%
%   SECD(F, PREF) does the same but uses the CHEBPREF object PREF when computing
%   the composition.
%
% See also ASECD, SEC.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Call the compose method:
F = compose(F, @secd, varargin{:});

end
