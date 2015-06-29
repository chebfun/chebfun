function F = cumprod(F)
%CUMPROD   Indefinite product integral.
%   CUMPROD(F) is the indefinite product integral of the CHEBFUN F, which is
%   defined as exp(cumsum(log(F))).
%
% See also PROD.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Loop over the columns:
for k = 1:numel(F)
    %[TODO]: Is this right? What about row CHEBFUNS?
    F(k) = exp(cumsum(log(F(k))));
end


end
