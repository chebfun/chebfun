function F = tan(F, varargin)
%TAN   Tangent of a CHEBFUN.
%   TAN(F) computes the tangent of the CHEBFUN F.
%
%   TAN(F, PREF) does the same but uses the CHEBFUNPREF object PREF when
%   computing the composition.
%
% See also ATAN, TAND.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% [TODO]:  Restore or change this once we have decided the proper behavior or
% isfinite() and defined that function.
% if ( ~isfinite(f) )
%     error('CHEBFUN:CHEBFUN:tan:inf',...
%         'SIN is not defined for functions which diverge to infinity');
% end

% Call the compose method:
F = compose(F, @tan, varargin{:});

end
