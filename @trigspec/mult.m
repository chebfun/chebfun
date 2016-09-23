function M = mult(disc, f)
%MULT    Multiplication operator for the TRIGSPEC class.
%   M = MULT(A, F) returns the multiplication operator that represents 
%   u(x) -> F(x)u(x), in the Fourier basis. 
% 
% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Obtaining some useful information:
n = disc.dimension;
dom = disc.domain;

% Convert to a TRIGTECH-based CHEBFUN. Use a non-adaptive construction: 
% we know that the TRIGTECH-based CHEBFUN will have a length smaller 
% than the CHEBTECH-based CHEBFUN.
f = chebfun(f, dom, 'periodic', length(f));
f = simplify(f);

% Call TRIGSPEC/MULTMAT.
M = trigspec.multmat(n, f);

end
