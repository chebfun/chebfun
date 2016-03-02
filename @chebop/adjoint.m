function Nstar = adjoint(N,u)
%ADJOINT   Compute the adjoint of a linear CHEBOP.
%   NSTAR = ADJOINT(N) returns the adjoint of a scalar CHEBOP that
%   has either periodic or endpoint constraints. If N is nonlinear then
%   N is first linearized around U = 0 and NSTAR is the adjoint of the 
%   linearization.
%
%   NSTAR = ADJOINT(N,U) linearizes around the function U then computes
%   the adjoint.
%
% See also LINOP/ADJOINT.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% set U to zero 
if ( nargin < 2 )
    u = 0*N.init;
end

% Linearize N around u.:
L = linearize(N, u);

% Get the value of the highest derivative:
n = L.diffOrder;

% Trivial case:
if ( n == 0 )
    Nstar = N;
    return
end

% call linop adjoint
[Lstar, op, bcOpL, bcOpR, bcOpM] = adjoint(L, getBCType(N));

%% 
% Construct chebop of the adjoint with pretty print bcs
Nstar = chebop(Lstar.domain);
Nstar.op = op;
Nstar.lbc = bcOpL;
Nstar.rbc = bcOpR;
Nstar.bc = bcOpM;

end

