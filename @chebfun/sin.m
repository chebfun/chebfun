function F = sin(F, varargin)
%SIN   Sine of a CHEBFUN.
%   SIN(F) computes the sine of the CHEBFUN F.
%
%   SIN(F, PREF) does the same but uses the CHEBPREF object PREF when computing
%   the composition.
%
% See also ASIN, SIND.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% [TODO]:  Restore or change this once we have decided the proper behavior or
% isfinite() and defined that function.
% if ( ~isfinite(f) )
%     error('CHEBFUN:sin:inf',...
%         'SIN is not defined for functions which diverge to infinity');
% end

% Call the compose method:
F = compose(F, @sin, varargin{:});

end
