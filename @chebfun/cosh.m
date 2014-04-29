function F = cosh(F, varargin)
%COSH   Hyperbolic cosine of a CHEBFUN.
%   COSH(F) computes the hyperbolic cosine of the CHEBFUN F.
%
%   COSH(F, PREF) does the same but uses the CHEBPREF object PREF when
%   computing the composition.
%
% See also COS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Call the compose method:
F = compose(F, @cosh, varargin{:});

end
