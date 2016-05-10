function F = sind(F, varargin)
%SIND   Sine of a CHEBFUN, result in degrees.
%   SIND(F) computes the sine (in degrees) of the CHEBFUN F.
%
%   SIND(F, PREF) does the same but uses the CHEBFUNPREF object PREF when computing
%   the composition.
%
% See also ASIND, SIN.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% [TODO]:  Restore or change this once we have decided the proper behavior or
% isfinite() and defined that function.
% if ( ~isfinite(f) )
%     error('CHEBFUN:CHEBFUN:sind:inf',...
%         'SIN is not defined for functions which diverge to infinity');
% end

% Call the compose method:
F = compose(F, @sind, varargin{:});

end
