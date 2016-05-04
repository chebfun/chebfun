function F = cscd(F, varargin)
%CSCD   Cosecant of a CHEBFUN, result in degrees.
%   CSCD(F) computes the cosecant (in degrees) of the CHEBFUN F.
%
%   CSCD(F, PREF) does the same but uses the CHEBFUNPREF object PREF when
%   computing the composition.
%
% See also ACSCD, CSC.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call the compose method:
F = compose(F, @cscd, varargin{:});

end
