function F = acsch(F, varargin)
%ACSCH   Inverse hyperbolic cosecant of a CHEBFUN.
%   ACSCH(F) computes the inverse hyperbolic cosecant of the CHEBFUN F.
%
%   ACSCH(F, PREF) does the same but uses the CHEBFUNPREF object PREF when
%   computing the composition.
%
% See also CSCH.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% Call the compose method:
F = compose(F, @acsch, varargin{:});

end
