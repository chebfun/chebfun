function F = csch(F, varargin)
%CSCH   Hyperbolic cosecant of a CHEBFUN.
%   CSCH(F) computes the hyperbolic cosecant of the CHEBFUN F.
%
%   CSCH(F, PREF) does the same but uses the CHEBFUNPREF object PREF when
%   computing the composition.
%
% See also ACSCH.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call the compose method:
F = compose(F, @csch, varargin{:});

end
