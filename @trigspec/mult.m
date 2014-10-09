function M = mult(disc, f)
%MULT    Multiplication operator for the TRIGSPEC class.
%   M = MULT(A, F) returns the multiplication operator that represents 
%   u(x) -> F(x)u(x), in the Fourier basis. 
% 
% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Obtaining some useful information:
n = disc.dimension;
d = disc.domain;

% Convert to a TRIGTECH-based CHEBFUN.
f = chebfun(f, d, 'periodic');

M = trigspec.multmat(n, f.funs{1});

end