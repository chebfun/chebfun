function F = sind(F, varargin)
%SIND   Sine of a CHEBFUN, result in degrees.
%   SIND(F) computes the sine (in degrees) of the CHEBFUN F.
%
%   SIND(F, PREF) does the same but uses the CHEBPREF object PREF when computing
%   the composition.
%
% See also ASIND, SIN.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% [TODO]:  Restore or change this once we have decided the proper behavior or
% isfinite() and defined that function.
% if ( ~isfinite(f) )
%     error('CHEBFUN:sin:inf',...
%         'SIN is not defined for functions which diverge to infinity');
% end

% Call the compose method:
F = compose(F, @sind, varargin{:});

end
