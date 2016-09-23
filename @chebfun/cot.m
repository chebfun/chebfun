function F = cot(F, varargin)
%COT   Cotangent of a CHEBFUN.
%   COT(F) computes the cotangent of the CHEBFUN F.
%
%   COT(F, PREF) does the same but uses the CHEBFUNPREF object PREF when
%   computing the composition.
%
% See also ACOT, COTD.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% [TODO]:  Restore or change this once we have decided the proper behavior or
% isfinite() and defined that function.
% if ( ~isfinite(f) )
%     error('CHEBFUN:CHEBFUN:cot:inf',...
%         'COT is not defined for functions which diverge to infinity');
% end

% Call the compose method:
F = compose(F, @cot, varargin{:});

end
