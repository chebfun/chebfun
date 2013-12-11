function F = acos(F, varargin)
%ACOS   Inverse cosine of a CHEBFUN.
%   ACOS(F) computes the inverse cosine of the CHEBFUN F.
%
%   ACOS(F, PREF) does the same but uses the CHEBPREF object PREF when
%   computing the composition.
%
% See also COS, ACOSD.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Call the compose method:
F = compose(F, @acos, varargin{:});

end
