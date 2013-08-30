function f = prod(f)
%PROD   Product integral.
%   CUMPROD(F) is the indefinite product integral of the CHEBFUN F, which is
%   defined as exp(sum(log(F))).

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f = exp(sum(log(f)));
%[TODO]: Is this right? What about row CHEBFUNS?

end