function F = tanh(F, varargin)
%TANH   Hyperbolic tangent of a CHEBFUN.
%   TANH(F) computes the hyperbolic tangent of the CHEBFUN F.
%
%   TANH(F, PREF) does the same but uses the CHEBFUNPREF object PREF when
%   computing the composition.
%
% See also ATANH.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call the compose method:
F = compose(F, @tanh, varargin{:});

end
