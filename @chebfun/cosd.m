function g = cosd(f, pref)
%COSD   Cosine of a chebfun, result in degrees.
%
% See also ACOSD, COS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfun.pref();
end

% [TODO]:  Restore or change this once we have decided the proper behavior or
% isfinite() and defined that function.
% if ( ~isfinite(f) )
%     error('CHEBFUN:cosd:inf',...
%         'COSD is not defined for functions which diverge to infinity');
% end

% Call the compose method:
g = compose(f, @cosd, pref);

end
