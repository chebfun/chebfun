function Nstar = adjoint(N, u)
%ADJOINT   Compute the adjoint of a CHEBOP.
%   NSTAR = ADJOINT(N) returns the adjoint of a CHEBOP that has either periodic
%   or endpoint boundary conditions. If N is nonlinear then N is first
%   linearized around U = 0 and NSTAR is the adjoint of the linearization.
%
%   NSTAR = ADJOINT(N, U) linearizes around the function U then computes
%   the adjoint.
%
% Examples:
%   L = chebop(-1,1); L.op = @(x,u) diff(u,2) + x*u;
%   L.lbc = 0; L.rbc = 'neumann', Ls = adjoint(L)
%
%   N = chebop(0,3); N.op = @(x,u) diff(u) + u^2;
%   u = chebfun('exp(x)',[0 3]);
%   N.lbc = 1, Ns = adjoint(N,u)
%
% See also LINOP/LINOPADJOINT.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Set U to zero:
if ( nargin < 2 )
    u = 0*N.init;
end

% Linearize N around u:
L = linearize(N, u);

% Get the value of the highest derivative:
n = max(max(L.diffOrder));

% Trivial case:
if ( n == 0 )
    Nstar = N;
    return
end

% Call linop adjoint:
[Lstar, op, bcOpL, bcOpR, bcOpM] = linopAdjoint(L, getBCType(N));

% Construct chebop of the adjoint with pretty print bcs:
Nstar = chebop(Lstar.domain);
Nstar.op = op;
Nstar.lbc = bcOpL;
Nstar.rbc = bcOpR;
Nstar.bc = bcOpM;

end