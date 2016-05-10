function F = csc(F, varargin)
%CSC   Cosecant of a CHEBFUN.
%   CSC(F) computes the cosecant of the CHEBFUN F.
%
%   CSC(F, PREF) does the same but uses the CHEBFUNPREF object PREF when
%   computing the composition.
%
% See also ACSC, CSCD.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call the compose method:
F = compose(F, @csc, varargin{:});

end
