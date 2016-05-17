function varargout = pde23t(varargin)
%PDE23TS   Solve PDEs using Chebfun.
%
%   UU = PDE23T(PDEFUN, TT, U0, BC) where PDEFUN is a handle to a function with
%   arguments u, t, x, and D, TT is a vector, U0 is a CHEBMATRIX, and BC is a
%   chebop boundary condition structure will solve the PDE dUdt = PDEFUN(UU, t,
%   x) with the initial condition U0 and boundary conditions BC over the time
%   interval TT.
%
%   This method is basically, a wrapper for @CHEBFUN/PDE23T(). See that file for
%   further details.
%
% See also CHEBFUN/PDE15S, PDESET.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% Convert a CHEBMATRIX input to a CHEBFUN:
varargin{3} = chebfun(varargin{3}); % Third input should be the only CHEBMATRIX.

[varargout{1:nargout}] = pde23t(varargin{:});

end
