function varargout = pde15s(varargin)
%PDE15S   Solve PDEs using the CHEBFUN system.
%   UU = PDE15s(PDEFUN, TT, U0, BC) where PDEFUN is a handle to a function with
%   arguments u, t, x, and D, TT is a vector, U0 is a CHEBMATRIX, and BC is a
%   chebop boundary condition structure will solve the PDE dUdt = PDEFUN(UU, t,
%   x) with the initial condition U0 and boundary conditions BC over the time
%   interval TT.
%
%   This method is basically, a wrapper for @CHEBFUN/PDES(). See that file for
%   further details.
%
% See also CHEBFUN/PDE15S

for j = 1:nargin
    if ( isa(varargin{j}, 'chebmatrix') )
        varargin{j} = chebfun(varargin{j});
    end
end

[varargout{1:nargout}] = pde15s(varargin{:});

end