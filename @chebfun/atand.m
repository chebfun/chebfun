function F = atand(F, varargin)
%ATAND   Inverse tangent of a CHEBFUN, result in degrees.
%   ATAND(F) computes the inverse tangent (in degrees) of the CHEBFUN F.
%
%   ATAND(F, PREF) does the same but uses the CHEBFUNPREF object PREF when
%   computing the composition.
%
% See also TAND, ATAN2D, ATAN.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% Call the compose method:
F = compose(F, @atand, varargin{:});

end
