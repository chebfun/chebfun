function Nstar = ctranspose(N)
%'   Compute the adjoint of a CHEBOP.
%   NSTAR = N' returns the adjoint of a CHEBOP that has either periodic or
%   endpoint boundary conditions. If N is nonlinear then N is first
%   linearized around U = 0 and NSTAR is the adjoint of the linearization.
%
%   ADJOINT(N) is called for the syntax N'.
%
% See also CHEBOP/ADJOINT.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% This is just a wrapper for adjoint().
Nstar = adjoint(N);

end
