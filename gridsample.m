function v = gridsample(f, n, dom, type)
%GRIDSAMPLE     Sample a function at gridpoints.
%   GRIDSAMPLE(F,N) returns F(CHEBPTS(N))
%   GRIDSAMPLE(F,N,DOM) returns F(CHEBPTS(N,DOM))
%   GRIDSAMPLE(F,N,DOM,'TRIG') returns F(TRIGPTS(N,DOM))

% This code has been written for Aurentz and Trefethen, 
% "Block operators and spectral discretizations", and is not
% yet fully compliant with standard Chebfun coding practices.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 2 )
    v = f(chebpts(n));
end

if ( nargin == 3 )
    v = f(chebpts(n,dom));
end

if ( nargin == 4 )
    v = f(trigpts(n,dom));
end

end
