function F = prod(F)
%PROD   Product integral.
%   PROD(F) is the definite product integral of the CHEBFUN F, which is defined
%   as exp(sum(log(F))).
%
% See also CUMPROD.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Loop over the columns:
for k = 1:numel(F)
    %[TODO]: Is this right? What about row CHEBFUNS?
    F(k) = exp(sum(log(F(k))));
end

end
