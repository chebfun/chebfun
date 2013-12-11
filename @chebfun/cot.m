function F = cot(F, pref)
%COT   Cotangent of a CHEBFUN.
%   COT(F) computes the cotangent of the CHEBFUN F.
%
%   COT(F, PREF) does the same but uses the CHEBPREF object PREF when
%   computing the composition.
%
% See also ACOT, COTD.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebpref();
end

% [TODO]:  Restore or change this once we have decided the proper behavior or
% isfinite() and defined that function.
% if ( ~isfinite(f) )
%     error('CHEBFUN:cot:inf',...
%         'COT is not defined for functions which diverge to infinity');
% end

% Call the compose method:
F = compose(F, @cot, pref);

end
