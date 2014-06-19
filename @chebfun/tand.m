function F = tand(F, varargin)
%TAND   Tangent of a CHEBFUN, result in degrees.
%   TAND(F) computes the tangent (in degrees) of the CHEBFUN F.
%
%   TAND(F, PREF) does the same but uses the CHEBFUNPREF object PREF when
%   computing the composition.
%
% See also ATAND, TAN.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call the compose method:
F = compose(F, @tand, varargin{:});

end
