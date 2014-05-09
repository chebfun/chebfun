function [L, res, isLinear] = linearize(N, u, x, flag)
%LINEARIZE    Linearize a CHEBOP.
%   L = LINEARIZE(N) returns a LINOP that corresponds to linearising the CHEBOP
%   N around the zero function on N.DOMAIN. The linop L will both include the
%   linearised differential equation, as well as linearised boundary conditions
%   (from N.LBC and N.RBC) and other constraints (from N.BC).
%
%   L = LINEARIZE(N, U) linearizes the CHEBOP N around the function U. Here, U
%   can either be a CHEBFUN or a CHEBMATRIX. IF U = [] then N is linearized
%   around the zero function, as above.
%
%   L = LINEARIZE(N, U, X) passes the independent variable, X, on N.DOMAIN.
%   If X = [] then LINEARIZE constructs the variable itself internally.
%
%   L = LINEARIZE(N, U, X, FLAG) is useful when we call LINOP(CHEBOP), i.e.,
%   converting a linear CHEBOP to a LINOP. If FLAG = 1, the method will stop
%   execution and return as soon as it encounters a nonlinear field in N. In
%   this case L is returned as an empty LINOP.
%
%   [L, RES] = LINEARIZE(N, ...) also returns RES; to the residual of the
%   differential equation part of N at the function it was linearized. In other
%   words, RES is the result of evaluating N at the zero function if no
%   additional function is passed to LINEARIZE(), or the function U if it is
%   passed. If N.OP is a scalar equation, RES is a CHEBFUN, otherwise it is a
%   CHEBMATRIX.
%
%   [L, RES, ISLINEAR] = LINEARIZE(N, ...) also returns the vector ISLINEAR,
%   with entries as follows:
%       ISLINEAR(1) = 1 if N.OP is linear, 0 otherwise.
%       ISLINEAR(2) = 1 if N.LBC is linear, 0 otherwise.
%       ISLINEAR(3) = 1 if N.RBC is linear, 0 otherwise.
%       ISLINEAR(4) = 1 if N.BC is linear, 0 otherwise.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% Start by assuming that N is linear.
isLinear = true(1, 4);

% Support single input argument for autonomous scalar problems:
if ( nargin(N) == 1 )
    N.op = @(x, u) N.op(u);
end

% Number of unknown variables N acts on. Subtract 1 from nargin(N.op), since the
% first argument is the independent variable x.
numVars = nargin(N.op) - 1;
dom = N.domain;

%% Construct a suitable function to linearize about:

% Construct the zero function on N.DOMAIN to linearize around if no U was
% passed.
if ( nargin < 2 || isempty(u) )
    % Initialise a zero CHEBFUN:
    zeroFun = chebfun(0, dom);
    % Wrap in a cell and call repmat() to get correct dimensions
    u = repmat({zeroFun}, numVars, 1);
end

% Construct the independent variable X if needed.
if ( nargin < 3 || isempty(x) )
    x = chebfun(@(x) x, dom);
end

% By default, set FLAG to 0.
if ( nargin < 4 || isempty(flag) )
    flag = 0;
end

% Convert the linerization variable to cell-array form:
if ( isa(u, 'chebmatrix') )
    u = u.blocks;
end
if ( isnumeric(u) ) 
    tmp = cell(numel(u), 1);
    for k = 1:numel(u)
    	tmp{k} = chebfun(u(k), N.domain);
    end
    u = tmp;
end
if ( ~iscell(u) )
    u = {u};
end

% Convert each element in the cell-array U to an ADCHEBFUN, and seed the
% derivative so it'll be of correct dimensions (i.e. correct block-size).
% Blocks corresponding to functions (i.e., CHEBFUNs) will be square, wheres the
% derivative blocks of scalars will by 1xINF.
isFun = ~cellfun(@isnumeric, u);
for k = 1:numVars
    u{k} = seed(adchebfun(u{k}, N.domain), k, isFun);
end

%% Evaluate N.op to get a linearisation of the differential equation:

% Evaluate N.op. The output will be the ADCHEBFUN NU. In case of systems, NU
% will be an array-valued ADCHEBFUN.
Nu = feval(N, x, u{:}); % N.op(x, u{:});

% Construct a LINOP L by vertically concatenating the derivatives stored in NU.
L = linop(vertcat(get(Nu, 'jacobian')));

% Construct the residual by vertically concatenating all functions stored in NU.
% RES will be a CHEBMATRIX.
res = vertcat(get(Nu, 'func'));

% Linearity information:
isLinear(1) = all(all(vertcat(get(Nu, 'linearity'))));

% If N is nonlinear, and we were looking to only test linearity, return.
if ( flag && ~all(isLinear) )
    L = linop();
    return
end

% Merge the domains of L obtained from evaluating the operator part above, with
% the domain of N, as we want to respect the breakpoints originally assigned
% to N (when its domain was defined)
L.domain = chebfun.mergeDomains(L.domain, dom);

%% Deal with parameterized problems:

% For problems with parameters, the system in L may not be square. This is OK if
% u0 contains doubles for the parameter entries. If not, we correct for this
% below by assuming the final few variables represent parameters.

[s1, s2] = size(L.blocks);
numParams = s2 - s1;
if ( all(isFun) && numParams > 0 )
    % We've found a paramterised problem, but weren't informed by u0. 
    
    % TODO: Do we really want to throw a warning?
%     % Throw a warning: 
%     if ( numParams == 1 )
%         warnStr = 'Assuming final variable is a parameter.';
%     else
%         warnStr = ['Assuming final ' int2str(numParams) ' variables are parameters.'];
%     end
%     warning('CHEBFUN:chebop:linearize:params', warnStr);
    
    % Reseed the final numParam variables as constants and linearize again:
    u = cellfun(@(b) b.func, u, 'UniformOutput', false);
    for k = 0:numParams-1
        u{end-k} = feval(u{end-k}, L.domain(1)); % Convert to a scalar.
    end
    [L, res, isLinear] = linearize(N, u, x, flag);
    return
end

%% Add BCs

% Initalise an empty LINOPCONSTRAINT.
BC = linopConstraint();

% Evaluate and linearise left boundary condition(s):
if ( ~isempty(N.lbc) )
    % Evaluate. The output, LBCU, will be an ADCHEBFUN.
    lbcU = N.lbc(u{:});
    
    % Ensure conditions were concatenated vertically, not horizontally
    lbcU = checkConcat(lbcU);
    
    % Loop through the components of LBCU.
    for k = 1:numel(lbcU)
        % Obtain the kth element of the ADCHEBFUN array.
        lbcUk = getElement(lbcU, k);
        % Evaluate the function at the left endpoint
        lbcUk = feval(lbcUk, dom(1));
        % Add the new condition to the LINOPCONSTRAINT BC.
        BC = append(BC, lbcUk.jacobian, lbcUk.func);
    end
    % Update linearity information.
    isLinear(2) = all(all(get(lbcU, 'linearity')));
end

% If N is nonlinear, and we were looking to only test linearity, return
if ( flag && ~all(isLinear) )
    L = linop();
    return
end

% Evaluate and linearise right boundary condition(s):
if ( ~isempty(N.rbc) )
    % Evaluate. The output, RBCU, will be an ADCHEBFUN.
    rbcU = N.rbc(u{:});
    
    % Ensure conditions were concatenated vertically, not horizontally
    rbcU = checkConcat(rbcU);
    
    % Loop through the components of RBCU.
    for k = 1:numel(rbcU)
        % Obtain the kth element of the ADCHEBFUN array.
        rbcUk = getElement(rbcU, k);
        % Evaluate the function at the right endpoint
        rbcUk = feval(rbcUk, dom(end));
        % Add the new condition to the LINOPCONSTRAINT BC.
        BC = append(BC, rbcUk.jacobian, rbcUk.func);
    end
    % Update linearity information.
    isLinear(3) = all(all(get(rbcU, 'linearity')));
end

% If N is nonlinear, and we were looking to only test linearity, return
if ( flag && ~all(isLinear) )
    L = linop();
    return
end

% Evaluate and linearise the remaining constraints:
if ( ~isempty(N.bc) )
    if ( strcmp(N.bc, 'periodic') )
        % Apply periodic boundary conditions:
       contConds = deriveContinuity(L, dom, true);
       contConds = contConds.continuity;
       BC = append(BC, contConds.functional, contConds.values);
       isLinear(4) = true;
    else
        % Evaluate. The output, BCU, will be an ADCHEBFUN.
        if ( nargin(N.bc) == 1 )
            bcU = N.bc(u{:});
        else
            bcU = N.bc(x, u{:});
        end

        % Ensure conditions were concatenated vertically, not horizontally
        bcU = checkConcat(bcU);

        % Gather all residuals of evaluating N.BC in one column vector.
        vals = cat(1, get(bcU, 'func'));
        % Loop through the conditions and append to the BC object.
        for k = 1:numel(bcU)
            BC = append(BC, get(bcU, 'jacobian', k), vals(k));
        end
        % Update linearity information.
        isLinear(4) = all(all(get(bcU, 'linearity')));
    end
end

% If N is nonlinear, and we were looking to only test linearity, return
if ( flag && ~all(isLinear) )
    L = linop();
    return
end

% if ( ~isempty(BC) && all(isFun) )
%     % Deal with parameterized problems:
%     BC.functional = diagonalise(BC.functional, isd);
% end

% Append all constraints to the LINOP returned.
L.constraint = BC;

end

function bc = checkConcat(bc)
% Ensure conditions were concatenated vertically, not horizontally.
if ( size(bc, 2) > 1 )
    warning('CHEBFUN:CHEBOP:linearize:bcConcat', ...
        ['CHEBOP conditions should be vertically concatenated with a '';'', ' ...
        'not horizontally with a '',''.\nHorizontal concatenation may ' ...
        'not be supported in future release.']);
    bc = bc.';
end
end