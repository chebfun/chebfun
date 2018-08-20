function W = mrdivide(V,c)
% MRDIVIDE / Right scalar divide for BALLFUNV object
%   MRDIVIDE(V,c) divides the BALLFUNV V by the scalar c

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

W = V*(1/c);
end
