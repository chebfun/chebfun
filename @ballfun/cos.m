function g = cos(f)
%COS   Cosine of a BALLFUN.
%   COS(F) computes the cosine of the BALLFUN F.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = f.coeffs;
G = ballfun.vals2coeffs(cos(ballfun.coeffs2vals(F)));
g = ballfun(G);
end
