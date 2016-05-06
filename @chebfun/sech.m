function F = sech(F, varargin)
%SECH   Hyperbolic secant of a CHEBFUN.
%   SECH(F) computes the hyperbolic secant of the CHEBFUN F.
%
%   SECH(F, PREF) does the same but uses the CHEBFUNPREF object PREF when computing
%   the composition.
%
% See also ASECH.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call the compose method:
F = compose(F, @sech, varargin{:});

end
