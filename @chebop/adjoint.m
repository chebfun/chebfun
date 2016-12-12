function Nstar = adjoint(N, u)
%ADJOINT   Compute the adjoint of a CHEBOP.
%   NSTAR = ADJOINT(N) returns the adjoint of a CHEBOP that has either periodic
%   or endpoint boundary conditions. If N is nonlinear then N is first
%   linearized around U = 0 and NSTAR is the adjoint of the linearization.
%
%   NSTAR = ADJOINT(N, U) linearizes around the function U then computes
%   the adjoint.
%
% See also LINOP/ADJOINT.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
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
[Lstar, op, bcOpL, bcOpR, bcOpM] = adjoint(L, getBCType(N));

% Construct chebop of the adjoint with pretty print bcs:
Nstar = chebop(Lstar.domain);
Nstar.op = op;
Nstar.lbc = bcOpL;
Nstar.rbc = bcOpR;
Nstar.bc = bcOpM;

end


function bcType = getBCType(N)
%GETBCTYPE   Detects the type of boundary conditions of a CHEBOP.
%   BCTYPE = GETBCTYPE(N) returns a string identifying the type of boundary
%   conditions used by a CHEBOP.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(N.bc) )
    bcType = 'bvp';
elseif ( strcmpi(N.bc, 'periodic') && isempty(N.lbc) && isempty(N.rbc) )
    bcType = 'periodic';
else
    bcType = 'general';
end

end
