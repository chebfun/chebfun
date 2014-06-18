function F = acot(F, varargin)
%ACOT   Inverse cotangent of a CHEBFUN.
%   ACOT(F) computes the inverse cotangent of the CHEBFUN F.
%
%   ACOT(F, PREF) does the same but uses the CHEBFUNPREF object PREF when
%   computing the composition.
%
% See also COT, ACOTD.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% Call the compose method:
F = compose(F, @acot, varargin{:});

end
