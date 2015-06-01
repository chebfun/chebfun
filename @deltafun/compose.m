function f = compose(f, op, g, data, pref)
%COMPOSE   Compose with a DELTAFUN is not allowed.
%
% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:DELTAFUN:compose:notSupported', ...
    'Composition with DELTAFUN is not supported.');

end
