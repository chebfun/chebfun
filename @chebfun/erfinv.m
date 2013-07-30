function g = erfinv(f, pref)
%ERFINV   Inverse error function of a chebfun.
%   X = ERFINV(Y) is the inverse error function for each element of Y. The
%   inverse error function satisfies Y = ERF(x), for -1 <= Y <= 1 and -inf <= X
%   <= inf.
%
% See also ERF, ERFC, ERFCX, ERFCINV.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Input must be real:
if ( ~isreal(f) )
    error('CHEBFUN:erfcinv:notreal', 'Input must be real.');
end

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfun.pref();
end

% Call the compose method:
g = compose(f, @erfinv, pref);

end