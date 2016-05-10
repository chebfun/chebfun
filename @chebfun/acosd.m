function F = acosd(F, varargin)
%ACOSD   Cosine of a CHEBFUN, result in degrees.
%   ACOSD(F) computes the cosine (in degrees) of the CHEBFUN F.
%
%   ACOSD(F, PREF) does the same but uses the CHEBFUNPREF object PREF when
%   computing the composition.
%
% See also ACOS, COS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% Call the compose method:
F = compose(F, @acosd, varargin{:});

end
