function g = cos(f)
% COS Cosine of a BALLFUN function
%   COS(f) is the cosine of the BALLFUN function f

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = f.coeffs;
G = ballfun.vals2coeffs(cos(ballfun.coeffs2vals(F)));
g = ballfun(G);
end
