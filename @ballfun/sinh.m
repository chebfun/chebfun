function g = sinh(f)
% SINH Hyperbolic sinus of a BALLFUN function
%   SINH(f) is the hyperbolic sinus of the BALLFUN function f

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = f.coeffs;
G = ballfun.vals2coeffs(sinh(ballfun.coeffs2vals(F)));
g = ballfun(G);
end
