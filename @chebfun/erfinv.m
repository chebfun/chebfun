function F = erfinv(F, pref)
%ERFINV   Inverse error function of a CHEBFUN.
%   X = ERFINV(Y) is the inverse error function of the CHEBFUN Y. The inverse
%   error function satisfies Y = ERF(x), for -1 <= Y <= 1 and -Inf <= X <= Inf.
%
% See also ERF, ERFC, ERFCX, ERFCINV.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Input must be real:
if ( ~isreal(F) )
    error('CHEBFUN:erfcinv:notreal', 'Input must be real.');
end

% Obtain preferences:
if ( nargin == 1 )
    pref = chebpref();
end

% Call the compose method:
F = compose(F, @erfinv, pref);

end
