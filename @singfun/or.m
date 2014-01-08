function f = and(f, g)
%|   SINGFUN logical OR is not supported.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

error('CHEBFUN:SINGFUN:or:notSupported',
    'Logical OR of two SINGFUNs is not supported.');

end
