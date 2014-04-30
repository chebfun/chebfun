function [L, res, isLinear] = linearize(N, u, x, flag)
%LINEARIZE    Linearize a CHEBOP.
%   L = LINEARIZE(N) returns a LINOP that corresponds to linearising the CHEBOP
%   N around the zero function on N.domain. The linop L will both include the
%   linearised differential equation, as well as linearised boundary conditions
%   (from N.LBC and N.RBC) and other constraints (from N.BC).
%
%   L = LINEARIZE(N, U) linearizes the CHEBOP N around the function U. Here, U
%   can either be a CHEBFUN or a CHEBMATRIX. IF U = [] then N is linearized
%   around the zero function, as above.
%
%   L = LINEARIZE(N, U, X) passes the independent variable, X, on N.domain.
%   If X = [] then LINEARIZE constructs the variable itself internally.
%
%   L = LINEARIZE(N, U, X, FLAG) is useful when we call LINOP(CHEBOP), i.e.,
%   converting a linear CHEBOP to a LINOP. If FLAG = 1, the method will stop
%   execution and return as soon as it encounters a nonlinear field in N. in
%   this case L is returned as an empty LINOP.
%
%   [L, RES] = LINEARIZE(N, ...) also returns the CHEBMATRIX RES, that
%   corresponds to the residual of the differential equation part of N at the
%   function it was linearized. In other words, RES is the result of evaluating
%   N at the zero function if no additioanl function is passed to LINEARIZE(),
%   or the function U if it is passed.
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

% To enable calling N.op, with a variable number of arguments, we actually want
% the function to be on a cell-array form.
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
if ( isa(u, 'chebfun') )
    u = {u};
end

% Convert each CHEBFUN object in the cell-array U to an ADCHEBFUN, and seed the
% derivative so it'll be of correct dimensions (i.e. correct block-size).
for k = 1:numVars
    u{k} = seed(adchebfun(u{k}, N.domain), k, numVars);
end

%% Evaluate N.op to get a linearisation of the differential equation:

% Evaluate N.op. The output will be the ADCHEBFUN NU. In case of systems, NU
% will be an array-valued ADCHEBFUN.
Nu = feval(N, x, u{:}); %N.op(x, u{:});

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
isd = isDiag(L);
L = diagonalise(L, isd);

%% Add BCs
% Initalise an empty LINOPCONSTRAINT.
BC = linopConstraint();

% Evaluate and linearise left boundary condition(s):
if ~( isempty(N.lbc) )
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
       contConds = contConds.continuity.functional;
       for k = 1:numel(contConds.blocks)
            BC = append(BC, contConds.blocks{k}, 0);
       end
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

if ( ~isempty(BC) )
    % Deal with parameterized problems:
    BC.functional = diagonalise(BC.functional, isd);
end

% Append all constraints to the LINOP returned.
L.constraint = BC;
end

function bc = checkConcat(bc)
% Ensure conditions were concatenated vertically, not horizontally
if ( size(bc, 2) > 1 )
    warning('CHEBFUN:CHEBOP:linearize:bcConcat', ...
        ['CHEBOP conditions should be vertically concatenated with a '';'', ' ...
        'not horizontally with a '',''.\nHorizontal concatenation may ' ...
        'not be supported in future release.']);
    bc = bc.';
end
end

function isd = isDiag(L)
% Start by assuming operators with zero difforder are diagonal operators:
isd = ~(L.diffOrder);
% Manually check these to confirm that they are really diagonal:
for k = find(isd).'
    tmp = chebmatrix(L.blocks(k));
    tmp = matrix(tmp, repmat(5, 1, numel(tmp.domain)-1));
    if ( norm(diag(diag(tmp))-tmp) > 1e-14 )
        isd(k) = false;
    end
end
end

function L = diagonalise(L, isd)
% Multiply the blocks in columns specified by ISD by the unitary CHEBFUN ONE.
blocks = L.blocks;
one = chebfun(1, L.domain);
for k = 1:size(blocks, 2)
    for j = 1:size(blocks,1);
        if ( all(isd(:,k)) )
            blocks{j,k} = blocks{j,k}*one;
        end
    end
end
L.blocks = blocks;
end