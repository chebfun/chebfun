function u = poisson(varargin)
%POISSON   Poisson solver in the unit ball with Dirichlet or Neumann
%          boundary conditions.
%   U = POISSON(F, BC, M, N, P) is the solution to the Poisson equation in
%   the unit ball with right-hand side F and Dirichlet boundary data
%   U(1,lambda,theta) = @(lambda,theta) BC(lambda, theta). It uses a
%   discretization size of M x N x P.
%
%   U = POISSON(F, BC, M) is the same as POISSON(F, BC, M, M, M).
%
%   U = POISSON(F, BC, M, N, P, 'neumann') is the solution to the Poisson 
%   equation with right-hand side F and Neuamnn boundary data U(1,lambda,theta) = 
%   @(lambda,theta) BC(lambda, theta). It uses a discretization size of M x N x P. 
%
%   U = POISSON(F, BC, M, 'neumann') is the same as POISSON(F, BC, M, M, M, 'neumann').
%
%   Also see HELMHOLTZ.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call Helmholtz command with zero frequency: 
u = helmholtz(varargin{1}, 0, varargin{2:end});
end