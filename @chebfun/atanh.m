function F = atanh(F, varargin)
%ATANH   Inverse hyperbolic tangent of a CHEBFUN.
%   ATANH(F) computes the inverse hyperbolic tangent of the CHEBFUN F.
%
%   ATANH(F, PREF) does the same but uses the CHEBFUNPREF object PREF when
%   computing the composition.
%
% See also TANH.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% Call the compose method:
F = compose(F, @atanh, varargin{:});

end
