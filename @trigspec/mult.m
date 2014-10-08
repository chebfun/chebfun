function M = mult(disc, f)
%MULT   Multiplication operator for the TRIGSPEC class.
%   M = MULT(A, F) returns the multiplication operator that represents 
%   u(x) -> F(x)u(x), in the Fourier basis. 
% 
% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Obtaining some useful information:
n = disc.dimension;
d = disc.domain;
f = restrict(f, d);
numIntervals = length(d) - 1;

% Find the diagonal blocks:
blocks = cell(numIntervals);
for k = 1:numIntervals
    blocks{k} = trigspec.multmat(n(k), f.funs{k});
end

% Assemble:
M = blkdiag(blocks{:});

end