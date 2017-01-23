function F = asech(F, varargin)
%ASECH   Inverse hyperbolic secant of a CHEBFUN.
%   ASECH(F) computes the inverse hyperbolic secant of the CHEBFUN F.
%
%   ASECH(F, PREF) does the same but uses the CHEBFUNPREF object PREF when
%   computing the composition.
%
% See also SECH.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% Call the compose method:
F = compose(F, @asech, varargin{:});

end
