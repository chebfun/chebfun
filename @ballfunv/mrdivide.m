function W = mrdivide(V,c)
%/    Right scalar divide for BALLFUNV object
%   V/c divides the BALLFUNV V by the scalar c.
%
% See also MLDIVIDE. 

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

W = V*(1/c);
end
