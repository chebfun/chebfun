function phi = phiFun(l)
%PHIFUN   Get a function handle to a phi-function.
%   PHI = PHIFUN(L) returns a function handle to the phi-function of index L.
%
% See also EXPINTEG/PHIEVAL, EXPINTEG/PSIFUN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Trivial case:
if ( l == 0 )
   phi = @(z) exp(z); 
   
% Compute them recursively using the recurrence formula 
% phi_{l}(z) = 1/z*(phi_{l-1}(z) - 1/(l-1)!).
else
   f = expinteg.phiFun(l-1);
   phi = @(z) (f(z) - 1/factorial(l-1))./z; 
end

end