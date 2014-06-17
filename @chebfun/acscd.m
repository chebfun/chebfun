function F = acscd(F, varargin)
%ACSCD   Inverse cosecant of a CHEBFUN, result in degrees.
%   ACSCD(F) computes the inverse cosecant (in degrees) of the CHEBFUN F.
%
%   ACSCD(F, PREF) does the same but uses the CHEBFUNPREF object PREF when
%   computing the composition.
%
% See also CSCD, ACSC.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% Call the compose method:
F = compose(F, @acscd, varargin{:});

end
