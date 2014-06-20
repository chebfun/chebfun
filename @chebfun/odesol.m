function varargout = odesol(sol, opt)
%ODESOL   Convert an ODE solution to CHEBFUN.
%   Y = ODESOL(SOL, OPT) converts the solution of an ODE initial-value or
%   boundary-value problem by standard MATLAB methods into a CHEBFUN
%   representation Y. SOL is the one-output form of any solver such as ODE45,
%   ODE15S, BVP5C, etc. OPT is the option structure used by the ODE solver. If
%   OPT is not passed, it is extracted from SOL.EXTDATA.OPTIONS if available.
%   The result is a piecewise CHEBFUN representing the solution.
%
%   [Y, T] = ODESOL(SOL, OPT) returns also the linear CHEBFUN T on the domain of
%   Y. Note that the order of outputs is the reverse of that from calls to
%   BVP4C(), ODE45(), etc.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%% Extract data from sol:
d = sol.x([1, end]);
vscale = max(abs(sol.y), [], 2); % Vertical scale (needed for RelTol)
numCols = size(sol.y, 1);

% Options:
if ( nargin < 2 ) 
    opt = [];
end
if ( isempty(opt) && isfield(sol, 'extdata') && ...
     isfield(sol.extdata, 'options') )
    % Take options from sol if none are given:
    opt = sol.extdata.options;
end

%% Find relative tolerances used in computations.
% Start with odeset default values:
relTol = 1e-6*ones(numCols, 1);         % Relative
absTol = 1e-6*ones(numCols, 1);         % Absolute
% Update if user used different tolerances:
if ( ~isempty(opt) )
    if ( ~isempty(opt.RelTol) )         % Relative tolerance given by user
        relTol = opt.RelTol*ones(numCols, 1);
    end
    if ( ~isempty(opt.AbsTol) )         % Absolute tolerance given by user
        if ( length(opt.AbsTol) == 1 )  % AbsTol might be vector or scalar
            absTol = opt.AbsTol*ones(numCols, 1);
        else
            absTol = opt.AbsTol;
        end
    end   
end
% Turn AbsTol into RelTol using scale:
relTol = max(relTol(:), absTol(:)./vscale(:));

%% Create a CHEBFUN object.
p = chebfunpref();
p.techPrefs.eps = max(relTol); % Use the same tolerance for each column.
p.splitting = true;            % use splitting, always, or there is no hope.
y = chebfun(@(x) deval(sol, x).', d, p);

% Parse outputs:
if ( nargout > 1 )
    t = chebfun('t', y.domain);
    varargout = {t, y};
else
    varargout = {y};
end

end

