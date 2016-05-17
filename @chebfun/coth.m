function F = coth(F, varargin)
%COTH   Hyperbolic cotangent of a CHEBFUN.
%   COTH(F) computes the hyperbolic cotangent of the CHEBFUN F.
%
%   COTH(F, PREF) does the same but uses the CHEBFUNPREF object PREF when
%   computing the composition.
%
% See also ACOTH.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call the compose method:
F = compose(F, @coth, varargin{:});

end
