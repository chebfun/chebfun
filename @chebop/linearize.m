function [L, res, isLinear, u] = linearize(N, u, x, linCheckFlag, paramReshapeFlag)
%LINEARIZE   Linearize a CHEBOP.
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
%   L = LINEARIZE(N, U, X, LINCHECKFLAG) is useful when we call LINOP(CHEBOP),
%   i.e., converting a linear CHEBOP to a LINOP. If LINCHECKFLAG = 1, the method
%   will stop execution and return as soon as it encounters a nonlinear field in
%   N. In this case L is returned as an empty LINOP. By default, LINCHECKFLAG =
%   0.
%
%   L = LINEARIZE(N, U, X, LINCHECKFLAG, PARAMRESHAPEFLAG) is useful when
%   determining whether parameters (as opposed to functions) appear in the
%   problem. If PARAMRESHAPEFLAG = 1, the code will try to cast any unknown
%   inputs which correspond to parameters to scalars, rather than CHEBFUNs, as
%   it linearizes. The Frechet derivatives corresponding to parameters will be
%   Inf x 1 CHEBFUNs, rather than Inf x Inf OPERATORBLOCKs. By default,
%   PARAMRESHAPEFLAG = 1. For generalized eigenvalue problems, it is useful to
%   pass PARAMRESHAPEFLAG = 0.
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
%
%   [L, RES, ISLINEAR, U] = LINEARIZE(N, ...) also returns CHEBMATRIX U that N
%   was linearized around. This is useful for parameter dependent problem, as
%   LINEARIZE() is where it is discovered that problems are parameter dependent,
%   so the CHEBMATRIX can be made to have to correct collection of CHEBFUN
%   objects and doubles, rather than just CHEBFUNs.
%
% See also LINOP.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% Start by assuming that N is linear.
isLinear = true(1, 4);

% The domain that the problem is specified on.
dom = N.domain;

%% Construct a suitable function to linearize about:

% Construct the zero function on N.DOMAIN to linearize around if no U was
% passed.
if ( nargin < 2 || isempty(u) )
    % Initialise a zero CHEBFUN:
    zeroFun = chebfun(0, dom);
    % Find out how many unknown variables N acts on.
    nVars = numVars(N);
    % Wrap in a cell and call repmat() to get correct dimensions
    u = repmat({zeroFun}, nVars, 1);
else
    if ( isa(u, 'chebmatrix') )
        nVars = size(u, 1);
    elseif ( isa(u, 'chebfun') )
        nVars = numColumns(u);
    else
        nVars = numel(u);
    end
end

% Construct the independent variable X if needed.
if ( nargin < 3 || isempty(x) )
    x = chebfun(@(x) x, dom);
end

% By default, set LINCHECKFLAG to 0.
if ( nargin < 4 || isempty(linCheckFlag) )
    linCheckFlag = 0;
end

% By default, set PARAMRESHAPEFLAG to 1.
if ( nargin < 5 || isempty(linCheckFlag) )
    paramReshapeFlag = 1;
end

% Convert the linearization variable to cell-array form:
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

% Support single input argument for autonomous scalar problems:
if ( nargin(N) == 1 )
    N.op = @(x, u) N.op(u);
end

% If nargin(N) == 2, but the dimension of the initial guess passed is greater
% than 1, we are working with the @(x,u) [diff(u{1}) + u{2}; ...] syntax. Need
% to make the code aware of this.
if ( nargin(N) == 2 && numel(u) > 1 )
    nVars = numel(u);
    cellArg = 1;
else
    cellArg = 0;
end

% Convert each element in the cell-array U to an ADCHEBFUN, and seed the
% derivative so it'll be of correct dimensions (i.e. correct block-size).
% Blocks corresponding to functions (i.e., CHEBFUNs) will be square, wheres the
% derivative blocks of scalars will by 1xINF.
isFun = ~cellfun(@isnumeric, u);
for k = 1:nVars
    u{k} = seed(adchebfun(u{k}, N.domain), k, isFun);
end

%% Evaluate N.op to get a linearisation of the differential equation:

% Evaluate N.op. The output will be the ADCHEBFUN NU. In case of systems, NU
% will be an array-valued ADCHEBFUN. We need different calling sequences
% depending on whether N has a cell-argument or not. We wrap the evaluation in a
% try-catch statement, since if we had an initial guess that leads to
% singularlity issues (such as N.init = 0 for N.op = @(u) diff(u,2) + sqrt(u)),
% we'd only throw meaningless error messages otherwise
try
    if ( cellArg )
        % No need to expand the cell U.
        Nu = feval(N, x, u);
    else
        % Need to expand the cell U.
        Nu = feval(N, x, u{:});
    end
catch ME
    if ( strcmp(ME.identifier, 'CHEBFUN:CHEBTECH:extrapolate:nansInfs') || ...
         strcmp(ME.identifier, ...
            'CHEBFUN:CHEBFUN:rdivide:columnRdivide:divisionByZeroChebfun') )
        error('CHEBFUN:CHEBOP:linearize:invalidInitialGuess', ...
            ['Failed to evaluate operator on the initial guess passed (or the ' ...
            'one constructed \nby CHEBOP). A potential cause might be ' ...
            'division by a zero CHEBFUN. Please supply\na valid initial ' ...
            'guess via the ''init'' field of the CHEBOP.'])
    else
        rethrow(ME)
    end
end
% Did the user specify the problem using old-school concatenation?
if ( size(Nu, 1) < size(Nu, 2) )
    warning('CHEBFUN:CHEBOP:linearize:vertcatOp', ...
        ['N.op should return a column vector.\n', ...
        'Row vectors are deprecated and may not be supported in future releases.'])
end

% Construct a LINOP L by vertically concatenating the derivatives stored in NU.
L = linop(vertcat(get(Nu, 'jacobian')));

% Construct the residual by vertically concatenating all functions stored in NU.
% RES will be a CHEBMATRIX.
res = vertcat(get(Nu, 'func'));

% Linearity information:
isLinear(1) = all(all(vertcat(get(Nu, 'linearity'))));

% If N is nonlinear, and we were looking to only test linearity, return.
if ( linCheckFlag && ~all(isLinear) )
    L = linop();
    return
end

% Merge the domains of L obtained from evaluating the operator part above, with
% the domain of N, as we want to respect the breakpoints originally assigned
% to N (when its domain was defined)
L.domain = domain.merge(L.domain, dom);

%% Deal with parameterized problems:

% For problems with parameters, the system in L may not be square. This is OK if
% u0 contains doubles for the parameter entries. If not, we correct for this
% below by assuming that any variable that does not have a diffOrder greater
% than 0 associated with it is a parameter, rather than a function. However, if
% no variable in the problem gets differentiated, we assume every variable is a
% function (unless we were told otherwise by an initial guess passed).
%
% The any(any(~L.isNotDiffOrInt)) checks if anything gets differentiated or
% integrated, the all(L.isNotDiffOrInt, 1) will then tell us what variables
% never were differentiated. The combination will only give us values of TRUE
% for variables that never get differentiated while other variables did.
isParam = any(any(~L.isNotDiffOrInt)) & all(L.isNotDiffOrInt, 1);

% If we have any parameters involved that are still thought to be functions, and
% we did not get a U passed in to linearize around, we reseed the corresponding
% variables.
if ( all(isFun) && any(isParam) && paramReshapeFlag )
    % We've found a parameterised problem, but weren't informed by u0.  Reseed
    % the final numParam variables as constants and linearize again:
    u = cellfun(@(b) b.func, u, 'UniformOutput', false);
    for k = find(isParam)
        u{k} = feval(u{k}, L.domain(1)); % Convert to a scalar.
    end
    [L, res, isLinear, u] = linearize(N, u, x, linCheckFlag);

    return
end

%% Add BCs.

% Initialise an empty LINOPCONSTRAINT.
BC = linopConstraint();

% Linearize left boundary condition:
if ( ~isempty(N.lbc) )
    [BC, isLinLeft] = linearizeLRbc(N.lbc, u, dom(1), BC, cellArg);
    isLinear(2) = isLinLeft;
end

% Linearize right boundary condition:
if ( ~isempty(N.rbc) )
    [BC, isLinRight] = linearizeLRbc(N.rbc, u, dom(end), BC, cellArg);
    isLinear(3) = isLinRight;
end

% Evaluate and linearise the remaining constraints. We need to treat the N.BC
% quite differently from N.LBC and N.RBC
if ( ~isempty(N.bc) )
    if ( strcmp(N.bc, 'periodic') )
       % Apply periodic boundary conditions:
       contConds = deriveContinuity(L, dom, true);
       contConds = contConds.continuity;
       BC = append(BC, contConds.functional, contConds.values);
       isLinear(4) = true;
    else
        % Before we evaluate the function, we need to ensure to contract the
        % derivative of any scalar ADchebfun. To see why, consider the boundary
        % condition for a problem on [0,1]:
        %   N.bc = @(x,v,p) v(0) - p;
        % where V is a CHEBFUN, but P is a scalar parameter. 
        %
        % When we reach this point of the code, the cell array U will contain 2
        % ADCHEBFUN objects. U{1} will be an ADCHEBFUN where U{1}.FUNC is a
        % CHEBFUN, and U{1}.JACOBIAN is a CHEBMATRIX that has the block types
        % [operatorBlock, chebfun]. U{2} will be an ADCHEBFUN where U{2}.FUNC is
        % a double, and U{2}.jacobian will be a CHEBMATRIX with the same block
        % types as U{1}.JACOBIAN.
        %
        % When we then evaluate the condition above, we first call FEVAL on V to
        % get the value V0 = V(0). This causes the dimensions of the jacobian of
        % V0 to collapse, so that V0.jacobian will be a CHEBMATRIX with the
        % block types [functionalBlock, double]. When we then do the subtraction
        % v(0) - p, the jacobian of P will be of incorrect dimensions,
        % reflecting that we haven't taken into account that we really are
        % evaluating a functional.
        %
        % To fix this, we hit the JACOBIAN field of any ADCHEBFUN currently in U
        % with an evaluation operator, which will cause the derivatives to
        % collapse to the correct dimensions, as required.
        if ( ~all(isFun) )
            % Create an evaluation operator for the left endpoint of the domain.
            % The evaluation point doesn't really matter, as all the derivatives
            % we're collapsing are either zeros or identity, so identical over
            % the whole domain
            E = functionalBlock.feval(dom(1), dom);
            for paramCounter = find(~isFun)
                u{paramCounter}.jacobian = E*u{paramCounter}.jacobian;
            end
        end
        
        % Evaluate. The output, BCU, will be an ADCHEBFUN.
        if ( nargin(N.bc) == 1 )
            bcU = N.bc(u{:});
        elseif ( cellArg )
            bcU = N.bc(x, u);
        else
            bcU = N.bc(x, u{:});
        end

        % Ensure conditions were concatenated vertically, not horizontally.
        bcU = checkConcat(bcU);

        % Gather all residuals of evaluating N.BC in one column vector.
        vals = cat(1, get(bcU, 'func'));
        
        % Loop through the conditions and append to the BC object.
        for k = 1:numel(bcU)
            J = get(bcU, 'jacobian', k);
            BC = append(BC, J , vals(k));
            jumps = get(bcU, 'jumpLocations', k);
            L = addGivenJumpAt(L,jumps);
        end
        
        % Update linearity information.
        isLinear(4) = all(all(get(bcU, 'linearity')));
        
        % Update domain:
        L.domain = domain.merge(L.domain, BC.functional.domain);
    end
    
end

% If N is nonlinear, and we were looking to only test linearity, return.
if ( linCheckFlag && ~all(isLinear) )
    L = linop();
    return
end

% Append all constraints to the LINOP returned.
L.constraint = BC;

% Cast the cell U back to a CHEBMATRIX, consisting of CHEBFUNs and scalars.
if ( nargout == 4)
    for k = 1:nVars
        u{k} = u{k}.func;
    end
    u = chebmatrix(u);
end

end

function [BC, isLinLR] = linearizeLRbc(op, u, evalPoint, BC, cellArg)
%LINEARIZELRBC   Linearize left and right boundary conditions.

% Evaluate N.lbc or N.rbc. The output will be the ADCHEBFUN LRBC. In case of
% systems, LRBC will be an array-valued ADCHEBFUN. We need different calling
% sequences depending on whether N has a cell-argument or not
if ( cellArg )
    % No need to expand the cell U
    lrBC = feval(op, u);
else
    % Need to expand the cell U
    lrBC = op(u{:});
end

% Ensure conditions were concatenated vertically, not horizontally
lrBC = checkConcat(lrBC);

% Loop through the components of LRBC.
for k = 1:numel(lrBC)
    % Obtain the kth element of the ADCHEBFUN array.
    lrBCk = getElement(lrBC, k);
    % Evaluate the function at the left endpoint
    lrBCk = feval(lrBCk, evalPoint);
    % Add the new condition to the LINOPCONSTRAINT BC.
    BC = append(BC, lrBCk.jacobian, lrBCk.func);
end

% Return linearity information
isLinLR = all(all(get(lrBC, 'linearity')));

end

function bc = checkConcat(bc)
% Ensure conditions were concatenated vertically, not horizontally.
if ( size(bc, 2) > 1 )
    warning('CHEBFUN:CHEBOP:linearize:bcConcat', ...
        ['CHEBOP conditions should be vertically concatenated with\n' ...
        'a '';'', not horizontally with a '',''. Horizontal concatenation may\n' ...
        'not be supported in future releases.']);
    bc = bc.';
end

end
