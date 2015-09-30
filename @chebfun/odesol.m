function varargout = odesol(sol, dom, opt)
%ODESOL   Convert an ODE solution to CHEBFUN.
% Y = ODESOL(SOL, DOM, OPT) converts the solution of an ODE initial-value or
% boundary-value problem by standard MATLAB methods into a CHEBFUN
% representation Y. The inputs to the method are:
%   SOL:   The one-output form of any solver such as ODE45, ODE15S, BVP5C, etc.
%          SOL can also be a cell-array of SOL structs, computed by resetting
%          the ODE solvers at breakpoints.
%   DOM:   The interval that the problem was solved on (may include
%          breakpoints).
%   OPT:   (Optional) the option structure used by the ODE solver. If OPT is not
%          passed, it is extracted from SOL.EXTDATA.OPTIONS if available.
%
% The output Y is a (potentially) piecewise CHEBFUN representing the solution.
%
%   [T, Y] = ODESOL(SOL, DOM, OPT) returns also the linear CHEBFUN T on the
%   domain of Y.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%% Extract data from sol:
% Compute vertical scale (needed for RelTol)
maxabs = @(sol) max(abs(sol.y), [], 2);
vscale = max(cell2mat(arrayfun(maxabs, sol, 'uniformOutput',false)), [], 2);

% Number of columns of the solution:
numCols = size(sol(1).y, 1);

% In the case where we have multiple entries in SOL, we construct a cell of
% function handles that we can evaluate to obtain a CHEBFUN. This is not needed
% in the single SOL case, and would actually cause an error later on if we have
% only one SOL for a piecewise domain (i.e. if restarting is turned off), so we
% treat the cases differently:
if ( length(sol) == 1 )
    devalFun = @(x) deval(sol, x).';
else
    dfun = @(sol) @(x) deval(sol, x).';
    % Obtain a cell of function handles that we can evaluate to obtain a
    % CHEBFUN:
    devalFun = arrayfun(dfun, sol, 'uniformOutput', false);
end

% Options:
if ( nargin < 3 ) 
    opt = [];
end
if ( isempty(opt) )
    % Take options from SOL if none are given. 
    if ( ~iscell(sol) && isfield(sol, 'extdata') && ...
        isfield(sol.extdata, 'options') )
        opt = sol.extdata.options;
    elseif ( iscell(sol) && isfield(sol{1}, 'extdata') && ...
            isfield(sol{1}.extdata, 'options') )
        % In multipiece case, take opts from the first piece. When solving IVPs
        % using CHEBOPs, OPTS should always get passed in from higher leves. 
        opt = sol{1}.extdata.options;
    end
end

% HappinessChecker
if ( ~isempty(opt) && isfield(opt, 'happinessCheck') )  
    checker = opt.happinessCheck;
else
    temp = cheboppref();
    checker = temp.happinessCheck;
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
idx = vscale == 0;
relTol(idx) = [];
absTol(idx) = [];
vscale(idx) = [];
relTol = max(relTol(:), absTol(:)./vscale(:));
if isempty(relTol)
    % We computed the zero solution, so relTol would be empty unless we give it
    % a value!
    relTol = eps;
end

%% Create a CHEBFUN object.
p = chebfunpref();
p.techPrefs.eps = max(relTol); % Use the same tolerance for each column.
p.techPrefs.happinessCheck = checker;
p.techPrefs.sampleTest = 0;
p.splitPrefs.splitMaxLength = 20000;

% Need to sort the domain D, since if we solve a final value problem, it will
% have been flipped.
[dom, idx] = sort(dom);

% If we flipped the domain, we also need to flip devalFun:
if ( idx(1) > idx(2) )
    devalFun = fliplr(devalFun);
end

% The output from DEVAL, based on what the ODE solvers return, is always
% vectorized, so we can turn the vectorcheck off.
y = chebfun(devalFun, dom, 'novectorcheck', p);

% Parse outputs:
if ( nargout > 1 )
    t = chebfun(@(t) t, y.domain);
    varargout = {t, y};
else
    varargout = {y};
end

end