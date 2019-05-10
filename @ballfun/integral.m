function g = integral(f, dim)
% INTEGRAL Definite integration of a BALLFUN.
%   INTEGRAL(F, DIM) where DIM is 1, 2 or 3 integrates only over r (radial direction), 
%   lambda (azimuthal direction) or theta (polar direction) respectively and 
%   and returns as its output a spherefun if DIM is 1 or a diskfun otherwise.
%
% See also SUM. 

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

g = sum(f, dim);
end