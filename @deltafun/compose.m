function f = compose(f, op, g, pref)
%COMPOSE   Compose with a DELTAFUN is not allowed.
%
% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error( 'CHEBFUN:DELTAFUN:COMPOSE', 'Composition with DELTAFUN is not supported.');

end
