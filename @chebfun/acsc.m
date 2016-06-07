function F = acsc(F, varargin)
%ACSC   Inverse cosecant of a CHEBFUN.
%   ACSC(F) computes the inverse cosecant of the CHEBFUN F.
%
%   ACSC(F, PREF) does the same but uses the CHEBFUNPREF object PREF when
%   computing the composition.
%
% See also CSC, ACSCD.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% Call the compose method:
F = compose(F, @acsc, varargin{:});

end
