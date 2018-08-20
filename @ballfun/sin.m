function g = sin(f)
% SIN Sinus of a BALLFUN function
%   SIN(f) is the sinus of the BALLFUN function f

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = f.coeffs;
G = ballfun.vals2coeffs(sin(ballfun.coeffs2vals(F)));
g = ballfun(G);
end
