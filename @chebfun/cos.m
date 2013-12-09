function F = cos(F, pref)
%COS   Cosine of a CHEBFUN.
%   COS(F) computes the cosine of the CHEBFUN F.
%
%   COS(F, PREF) does the same but uses the CHEBPREF object PREF when
%   computing the composition.
%
% See also ACOS, COSD.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebpref();
end

% [TODO]:  Restore or change this once we have decided the proper behavior or
% isfinite() and defined that function.
%if ( ~isfinite(f) )
%    error('CHEBFUN:cos:inf',...
%        'COS is not defined for functions which diverge to infinity');
%end

% Loop over the columns of F:
for k = 1:numel(F)
    % Call the compose method:
    F(k) = compose(F(k), @cos, pref);
end

end
