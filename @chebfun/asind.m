function F = asind(F, varargin)
%ASIND   Inverse sine of a CHEBFUN, result in degrees.
%   ASIND(F) computes the inverse sine (in degrees) of the CHEBFUN F.
%
%   ASIND(F, PREF) does the same but uses the CHEBFUNPREF object PREF when
%   computing the composition.
%
% See also SIND, ASIN.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% Call the compose method:
F = compose(F, @asind, varargin{:});

end
