function f = cumprod(f)
%CUMPROD   Indefinite product integral.
%   CUMPROD(F) is the indefinite product integral of the CHEBFUN F, which is
%   defined as exp(cumsum(log(F))).

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f = exp(cumsum(log(f)));
%[TODO]: Is this right? What about row CHEBFUNS?

end