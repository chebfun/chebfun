function F = sinh(F, varargin)
%SINH   Hyperbolic sine of a CHEBFUN.
%   SINH(F) computes the hyperbolic sine of the CHEBFUN F.
%
%   SINH(F, PREF) does the same but uses the CHEBFUNPREF object PREF when computing
%   the composition.
%
% See also ASINH.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call the compose method:
F = compose(F, @sinh, varargin{:});

end
