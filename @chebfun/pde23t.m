function varargout = pde23t(varargin)
%PDE23T   Solve PDEs using Chebfun.
%
%   PDE23T() has identical functionality to PDE15S(), but often performs better
%   than PDE15S() for non-diffusive problems. See the PDE15() help text for
%   further information.
%
% See also PDESET, PDE15S, ODE23T, PDESOLVE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% Loop over the inputs and look for a PDESET struct:
optsPassed = false;                 % Set flag if we find optsd were passed.
for k = 1:nargin
    if ( isstruct(varargin{k}) && isfield(varargin{k}, 'ODESolver') )
        % Ammend the ODESolver field to be @ode23t:
        varargin{k}.ODESolver = @ode23t;
        optsPassed = true;          % We know now that options were passed.
        break                       % No need to continue the loop.
    end
end

% If no PDESET options were passed, create a PDESET struct:
if ( ~optsPassed )
    opts = pdeset('ODESolver', @ode23t);
    varargin{end+1} = opts;
end

% Call PDESOLVE() with option to use ODE23T():
[varargout{1:nargout}] = pdeSolve(varargin{:});

end
