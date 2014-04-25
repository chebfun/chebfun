function [L, res, isLinear] = linearize(N, x, u, flag)
% LINEARIZE     Linearize a CHEBOP.
%
% L = LINEARIZE(N) returns a LINOP that corresponds to linearising the CHEBOP N
%   around the zero function on N.domain. The linop L will both include the
%   linearised differential equation, as well as linearised boundary conditions
%   (from N.LBC and N.RBC) and other constraints (from N.BC).
%
% L = LINEARIZE(N, X, U) linearizes the CHEBOP N around the function U. Here, X
%   denotes the independent variable on N.domain. Note that calling
%   LINEARIZE(N, [], U) is also supported, i.e., the method will create the
%   variable X. Here, U can either be a CHEBFUN or a CHEBMATRIX.
%
% L = LINEARIZE(N, X, U, FLAG) is useful when we call LINOP(CHEBOP), i.e.
%   converting a linear CHEBOP to a LINOP. If FLAG = 1, the method will stop
%   execution and return as soon as it encounters a nonlinear field in N.
%
% [L, RES] = LINEARIZE(N, ...) also returns the CHEBMATRIX RES, that corresponds
%   to the residual of the differential equation part of N at the function it
%   was linearized. In other words, RES is the result of evaluating N at the
%   zero function if no function is passed to LINEARIZE(), or the function U if
%   it is passed.
%
% [L, RES, ISLINEAR] = LINEARIZE(N, ...) also returns the vector ISLINEAR, with
%   entries as follows:
%       ISLINEAR(1) = 1 if N.OP is linear, 0 otherwise.
%       ISLINEAR(2) = 1 if N.LBC is linear, 0 otherwise.
%       ISLINEAR(3) = 1 if N.RBC is linear, 0 otherwise.
%       ISLINEAR(4) = 1 if N.BC is linear, 0 otherwise.

% TODO: Flag should determine whether we just want to do a linearity check, i.e.
% stop if we're trying to convert a chebop to a linop.

% Start by assuming that N is linear.
isLinear = true(1, 4);
% Number of unknown variables N acts on. Subtract 1 from nargin(N.op), since the
% first argument is the independent variable x.
numVars = nargin(N.op) - 1;

% Construct the independent variable X if needed.
if ( nargin < 2 || isempty(x) )
    x = chebfun(@(x) x, N.domain);
end

% Construct the zero function on N.DOMAIN to linearize around if no U was
% passed.
if ( nargin < 3 || isempty(u) )
    % Initialise a zero CHEBFUN:
    zeroFun = chebfun(0, N.domain);
    % Wrap in a cell and call repmat() to get correct dimensions
    u = repmat({zeroFun}, numVars, 1);
end

% By default, set FLAG to 0.
if ( nargin < 4 || isempty(flag) )
    flag = 0;
end

% To enable calling N.op, with a variable number of arguments, we actually want
% the function to be on a cell-array form.
if ( isa(u, 'chebmatrix') )
    u = u.blocks;
end
if ( isa(u, 'chebfun') )
    u = {u};
end

% Convert each CHEBFUN object in the cell-array U to an ADCHEBFUN, and seed the
% derivative so it'll be of correct dimensions (i.e. correct block-size).
for k = 1:numVars
    u{k} = seed(adchebfun(u{k}), k, numVars);
end


%% Evaluate N.op to get a linearisation of the differential equation:

% Evaluate N.op. The output will be the ADCHEBFUN NU. In case of systems, NU
% will be an array-valued ADCHEBFUN.
Nu = N.op(x, u{:});

% Construct the LINOP L by vertically concatenating the derivatives stored in
% NU.
L = linop(vertcat(get(Nu, 'jacobian')));
% Construct the residual by vertically concatenating all functions stored in NU.
% RES will be a CHEBMATRIX.
res = vertcat(get(Nu, 'func'));
% Linearity information
isLinear(1) = all(all(vertcat(get(Nu, 'linearity'))));

% If N is nonlinear, and we were looking to only test linearity, return.
if ( flag && ~all(isLinear) )
    return
end

% Merge the domains of L obtained from evaluating the operator part above, with
% the domain of N, as we want to respect the breakpoints originally assigned
% to N (when its domain was defined)
L.domain = chebfun.mergeDomains(L.domain, N.domain);

%% Add BCs
% Initalise an empty LINOPCONSTRAINT.
BC = linopConstraint();

% Evaluate and linearise left boundary condition(s):
if ~( isempty(N.lbc) )
    % Evaluate. The output, LBCU, will be an ADCHEBFUN.
    lbcU = N.lbc(u{:});
    
    % Ensure conditions were concatenated vertically, not horizontally
    checkConcat(lbcU);
    
    % Loop through the components of LBCU.
    for k = 1:numel(lbcU)
        % Obtain the kth element of the ADCHEBFUN array.
        lbcUk = getElement(lbcU, k);
        % Evaluate the function at the left endpoint
        lbcUk = feval(lbcUk, N.domain(1));
        % Add the new condition to the LINOPCONSTRAINT BC.
        BC = append(BC, lbcUk.jacobian, lbcUk.func);
    end
    % Update linearity information.
    isLinear(2) = all(all(get(lbcU, 'linearity')));
end

% If N is nonlinear, and we were looking to only test linearity, return
if ( flag && ~all(isLinear) )
    return
end

% Evaluate and linearise right boundary condition(s):
if ( ~isempty(N.rbc) )
    % Evaluate. The output, RBCU, will be an ADCHEBFUN.
    rbcU = N.rbc(u{:});
    
    % Ensure conditions were concatenated vertically, not horizontally
    checkConcat(rbcU);
    
    % Loop through the components of RBCU.
    for k = 1:numel(rbcU)
        % Obtain the kth element of the ADCHEBFUN array.
        rbcUk = getElement(rbcU, k);
        % Evaluate the function at the right endpoint
        rbcUk = feval(rbcUk, N.domain(end));
        % Add the new condition to the LINOPCONSTRAINT BC.
        BC = append(BC, rbcUk.jacobian, rbcUk.func);
    end
    % Update linearity information.
    isLinear(3) = all(all(get(rbcU, 'linearity')));
end

% If N is nonlinear, and we were looking to only test linearity, return
if ( flag && ~all(isLinear) )
    return
end

% Evaluate and linearise the remaining constraints:
if ( ~isempty(N.bc) )
    % Evaluate. The output, BCU, will be an ADCHEBFUN.
    bcU = N.bc(x, u{:});
    
    % Ensure conditions were concatenated vertically, not horizontally
    checkConcat(bcU);
    
    % Gather all residuals of evaluating N.BC in one column vector.
    vals = cat(1, get(bcU, 'func'));
    % Loop through the conditions and append to the BC object.
    for k = 1:numel(bcU)
        BC = append(BC, get(bcU, 'jacobian', k), vals(k));
    end
    % Update linearity information.
    isLinear(4) = all(all(get(bcU, 'linearity')));
end
% Append all constraints to the LINOP returned.
L.constraint = BC;
end

function checkConcat(bc)
% Ensure conditions were concatenated vertically, not horizontally
if ( size(bc, 2) > 1 )
    error('CHEBFUN:CHEBOP:linearize:bcConcat', ...
        ['Chebop conditions must be vertically concatenated with a '';'',' ...
        '\nnot horizontally with a '',''.']);
end
end