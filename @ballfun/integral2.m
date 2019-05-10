function g = integral2(f, dims)
% INTEGRAL2 Definite integration of a BALLFUN in two variables.
%   INTEGRAL2(F, DIMS) integrates F over two of the variables r, lambda or theta 
%   where DIMS is a row vector containing two of the three indices 1, 2 or 3.
%   The output is a 1D CHEBFUN in the remaining variable.
%
%   G = INTEGRAL2(F) is the same as INTEGRAL2(F, [2, 3]).
%
% See also SUM2. 

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

g = sum2(f, dims);
end