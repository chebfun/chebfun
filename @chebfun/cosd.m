function g = cosd(f, pref)
%COSD   Cosine of a CHEBFUN, result in degrees.
%   COSD(F) computes the cosine (in degrees) of the CHEBFUN F.
%
%   COSD(F, PREF) does the same but uses the CHEBPREF object PREF when
%   computing the composition.
%
% See also ACOSD, COS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebpref();
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
