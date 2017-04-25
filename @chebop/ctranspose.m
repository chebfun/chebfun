function Nstar = ctranspose(N)
%'   Compute the adjoint of a CHEBOP.
%   NSTAR = N' returns the adjoint of a CHEBOP that has either periodic or
%   endpoint boundary conditions. If N is nonlinear then N is first
%   linearized around U = 0 and NSTAR is the adjoint of the linearization.
%
%   N' is calculated by a call to ADJOINT(N).  To linearize about
%   a nonzero function U, use ADJOINT(N, U).
%
%  Example:
%   L = chebop(-1,1); L.op = @(x,u) diff(u,2) + diff(u,1) + u;
%   L.lbc = 0; L.rbc = 1, Ls = L'
%
% See also CHEBOP/ADJOINT.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% This is just a wrapper for adjoint().
Nstar = adjoint(N);

end
