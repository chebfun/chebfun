function F = outer(disc, f, g)
%OUTER   Outer product operator in VALSDISCRETIZATION.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Obtain the collocation points and associated weights.
[x, w] = functionPoints(disc);

% Evaluate F and G at the grid:
fx = feval(f, x);
gx = feval(g, x);

% Repmat the weights W if needed (in case of array valued CHEBFUNs)
w = repmat(w, size(fx, 2), 1);
% Take the outer product to obtain the matrix we desire.
F = fx * (w .* gx);

end
