function F = acotd(F, varargin)
%ACOTD   Inverse cotangent of a CHEBFUN, result in degrees.
%   ACOTD(F) computes the inverse cotangent (in degrees) of the CHEBFUN F.
%
%   ACOTD(F, PREF) does the same but uses the CHEBFUNPREF object PREF when
%   computing the composition.
%
% See also COTD, ACOT.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% Call the compose method:
F = compose(F, @acotd, varargin{:});

end
