function Nstar = adjoint(N)
%ADJOINT   Compute the adjoint of a linear CHEBOP.
%   ADJOINT(N), where N is a CHEBOP, returns the adjoint CHEBOP of N.
%
%   NSTAR = ADJOINT(N) also returns the adjoint of a scalar CHEBOP that
%   has either periodic or endpoint constraints.
%
% See also LINOP/ADJOINT.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Tell CHEBOP/LINEARIZE to stop if it detects nonlinearity:
linCheck = true; 

% Linearize, thereby obtaining linearity information, a LINOP, and an input of
% the correct dimensions to pass to N:
[L, ~, isLinear, ~] = linearize(N, N.init, [], linCheck);

% We need the entire operator (including BCs) to be linear:
assert(all(isLinear), 'CHEBFUN:CHEBOP:adjoint:nonlinear', ...
    ['The input operator appears to be nonlinear.\n', ...
    'ADJOINT supports only linear CHEBOP instances.']);

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

