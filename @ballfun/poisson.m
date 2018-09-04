function u = poisson(f, g, m, n, p)
%POISSON   Poisson solver with Dirichlet boundary conditions.
%   POISSON(F, G, m, n, p) is the solution to the Poisson
%   equation with right-hand side F and Dirichlet boundary
%   data given by g(lambda, theta).
%
% Also see HELMHOLTZ, HELMHOLTZ_NEUMANN.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call Helmholtz command with zero frequency: 
u = helmholtz(f, 0, g, m, n, p);

end