function F = asinh(F, varargin)
%ASINH   Inverse hyperbolic sine of a CHEBFUN.
%   ASINH(F) computes the inverse hyperbolic sine of the CHEBFUN F.
%
%   ASINH(F, PREF) does the same but uses the CHEBFUNPREF object PREF when
%   computing the composition.
%
% See also SINH.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% Call the compose method:
F = compose(F, @asinh, varargin{:});

end
