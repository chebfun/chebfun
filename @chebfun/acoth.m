function F = acoth(F, varargin)
%ACOTH   Inverse hyperbolic cotangent of a CHEBFUN.
%   ACOTH(F) computes the inverse hyperbolic cotangent of the CHEBFUN F.
%
%   ACOTH(F, PREF) does the same but uses the CHEBFUNPREF object PREF when
%   computing the composition.
%
% See also COTH.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% Call the compose method:
F = compose(F, @acoth, varargin{:});

end
