function g = sqrt(f)
%SQRT   Square root of a BALLFUN.
%   SQRT(F) is the square root of the BALLFUN F. 

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = f.coeffs;
G = ballfun.vals2coeffs(sqrt(ballfun.coeffs2vals(F)));
g = ballfun(G);
end
