function F = gamma(F, varargin)
%GAMMA   Gamma function of a CHEBFUN.
%   GAMMA(F) computes the composition of the gamma function with F.
%
%   GAMMA(F, PREF) does the same but uses the CHEBFUNPREF object PREF when
%   computing the composition.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% [TODO]:  Restore or change this once we have decided the proper behavior or
% isfinite() and defined that function.
%if ( ~isfinite(f) )
%    error('CHEBFUN:CHEBFUN:cos:inf',...
%        'COS is not defined for functions which diverge to infinity');
%end

% Call the compose method:
F = compose(F, @gamma, varargin{:});

end
