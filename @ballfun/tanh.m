function g = tanh(f)
% TANH Hyperbolic tangent of a BALLFUN function
%   TANH(f) is the hyperbolic tangent of the BALLFUN function f

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = f.coeffs;
G = ballfun.vals2coeffs(tanh(ballfun.coeffs2vals(F)));
g = ballfun(G);
end
