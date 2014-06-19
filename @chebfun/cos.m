function F = cos(F, varargin)
%COS   Cosine of a CHEBFUN.
%   COS(F) computes the cosine of the CHEBFUN F.
%
%   COS(F, PREF) does the same but uses the CHEBFUNPREF object PREF when
%   computing the composition.
%
% See also ACOS, COSD.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% [TODO]:  Restore or change this once we have decided the proper behavior or
% isfinite() and defined that function.
%if ( ~isfinite(f) )
%    error('CHEBFUN:CHEBFUN:cos:inf',...
%        'COS is not defined for functions which diverge to infinity');
%end

% Call the compose method:
F = compose(F, @cos, varargin{:});

end
