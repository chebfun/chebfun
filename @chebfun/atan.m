function F = atan(F, varargin)
%ATAN   Inverse tangent of a CHEBFUN.
%   ATAN(F) computes the inverse tangent of the CHEBFUN F.
%
%   ATAN(F, PREF) does the same but uses the CHEBFUNPREF object PREF when
%   computing the composition.
%
% See also TAN, ATAND.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% Call the compose method:
F = compose(F, @atan, varargin{:});

end
