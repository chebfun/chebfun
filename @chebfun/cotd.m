function F = cotd(F, varargin)
%COSD   Cotangent of a CHEBFUN, result in degrees.
%   COSD(F) computes the cotangent (in degrees) of the CHEBFUN F.
%
%   COSD(F, PREF) does the same but uses the CHEBFUNPREF object PREF when
%   computing the composition.
%
% See also ACOTD, COT.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% [TODO]:  Restore or change this once we have decided the proper behavior or
% isfinite() and defined that function.
% if ( ~isfinite(f) )
%     error('CHEBFUN:CHEBFUN:cotd:inf',...
%         'COTD is not defined for functions which diverge to infinity');
% end

% Call the compose method:
F = compose(F, @cotd, varargin{:});

end
