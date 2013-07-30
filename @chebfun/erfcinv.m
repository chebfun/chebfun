function g = erfcinv(f, pref)
%ERFCINV   Inverse complementary error function of a chebfun.
%   X = erfcinv(Y) is the inverse of the complementary error function for each
%   element of Y. It satisfies Y = ERFC(X) for 2 >= Y >= 0 and -Inf <= X <=
%   Inf.
%
% See also ERF, ERFC, ERFCX, ERFINV.

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
g = compose(f, @erfcinv, pref);

end