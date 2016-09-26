function F = asin(F, varargin)
%ASIN   Inverse sine of a CHEBFUN.
%   ASIN(F) computes the inverse sine of the CHEBFUN F.
%
%   ASIN(F, PREF) does the same but uses the CHEBFUNPREF object PREF when
%   computing the composition.
%
% See also SIN, ASIND.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% Call the compose method:
F = compose(F, @asin, varargin{:});

end
