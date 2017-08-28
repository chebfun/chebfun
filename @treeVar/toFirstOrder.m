function [funOut, indexStart, problemDom, coeffs, totalDiffOrders] = ...
    toFirstOrder(funIn, rhs, domain, numArgs, cellArg)
%TOFIRSTORDER   Convert higher order anonymous functions to first order systems.
%   Calling sequence:
%      [FUNOUT, INDEXSTART, PROBLEMDOM, COEFFS, TOTALDIFFORDERS] = ...
%          TOFIRSTORDER(FUNIN, RHS, DOMAIN)
%   where the inputs are:
%      FUNIN:  An anonymous function which describes an ODE, usually including
%              higher order derivatives, e.g. @(x,u) diff(u, 2) + sin(u).
%              Usually, the anonymous function comes from the OP field of a
%              CHEBOP.
%      RHS:    The right hand side of the differential equation being solved
%              with a CHEBOP, that is, the right argument of a CHEBOP
%              backslash.
%      DOMAIN: The domain of the problem we're trying to solve.
%   and the outputs are:
%      FUNOUT:     An anonymous function that is the first order reformulation
%                  of FUNIN, which can be passed to the MATLAB solvers.
%      INDEXSTART: A vector that denotes at which index we should start
%                  indexing each variable from so that they're in the correct
%                  order for the MATLAB ODE solvers.
%      PROBLEMDOM: The domain on which the problem can be specified, which may
%                  include breakpoints originally not included in DOMAIN.
%      COEFFS:     A cell array of the coefficients the multiply the highest
%                  order derivatives in problem. For example, for the problem
%                  @(x,u) 5*diff(u, 2) + u, we will have COEFFS{1} = 5.
%      TOTALDIFFORDERS:
%                  A vector that contains the maximum diffOrder applying to
%                  each variable in the problem.
%
%   In addition, one can pass the following arguments to the function, which are
%   only required if working in CHEBMATRIX syntax, for example, 
%   @(x,u)[diff(u{1}) + u{2};...]:
%       NVARS:   The number of unknown variables in the ODE.
%       CELLARG: Specify as TRUE if FUNIN is specified in CHEBMATRIX syntax.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Independent variable on the domain
t = chebfun(@(t) t, domain);

% See if we need to determine the number of variables involved, or if it was
% passed as an argument.
if ( nargin < 4)
    % We always have at least one argument to FUNIN, even in the scalar case. In
    % the system case, the first argument to funIn must be the independent time
    % variable. So the number of treeVar arguments needed is one less:
    numArgs = max(1, nargin(funIn) - 1);
    % If working in CHEBMATRIX syntax, we need all possible input variables to
    % be passed in as arguments:
    cellArg = false;
end
args = cell(numArgs, 1);
argsVec = zeros(1, numArgs);

% Populate the args cell
for argCount = 1:numArgs
    argsVec(argCount) = 1;
    args{argCount} = treeVar(argsVec, domain);
    % Reset the index vector
    argsVec = 0*argsVec;
end

% Evaluate FUNIN with the TREEVAR arguments. Need different calling sequences
% depending on whether we're working with CHEBMATRIX syntax and cells or not:
if ( cellArg )
    fevalResult = funIn(t, args);
else
    % If FUNIN only takes one argument, we just give it the TREEVAR argument. 
    % Otherwise, the first input will be the independent variable on the domain:
    if ( nargin(funIn) == 1 )
        fevalResult = funIn(args{:});
    else
        fevalResult = funIn(t, args{:});
    end
end

% If we got passed the problem as N\0, i.e. a RHS that did not match the
% dimensions, we need to repmat it so that the RHS variable has the correct
% dimensions below:
if ( isnumeric(rhs) && length(rhs) == 1 )
    rhs = repmat(rhs, size(fevalResult));
end

% Ensure RHS is a CHEBMATRIX
if ( ~isa(rhs, 'chebmatrix') )
    rhs = chebmatrix(rhs, domain);
end


% Initialize cells to store the infix forms of the expressions, the coefficients
% multiplying the highest order derivatives and any variables that appear in the
% anonymous functions:
systemInfix = cell(length(fevalResult), 1);
coeffs = systemInfix;
varArrays = systemInfix;

% We want to return all potential breakpoints caused by piecewise coefficients
% in the problem. So loop through fevalResult, and take the union of the
% breakpoints:
problemDom = fevalResult(1).domain;
for resCounter = 2:length(fevalResult)
    problemDom = union(problemDom, fevalResult(resCounter).domain);
end

% Ensure we also include potential breakpoints from the RHS:
problemDom = union(problemDom, rhs.domain);

% First look at all diffOrders to ensure we start with the correct indices.
% INDEXSTART denotes at which index we should start indexing each variable from.
% E.g., in the coupled system
%   [v'' + w; v + w'']
% we will have v = u(1), v' = u(2), w = u(3), w' = u(4), so INDEXSTART = [1, 3],
% since v and its derivatives starts getting counted at 1, and w and its
% derivatives start getting counted at 3.
%
% The vector INDEXSTARTDER is similar, but here, we also assign the highest
% order derivative of each variable it's own index. This is so that later on, we
% can correctly evaluate the coefficient multiplying the highest order
% derivative in problem. Thus, in the example above, we'd have
%   v = u(1), v' = u(2), v'' = u(3), w = u(4), w' = u(5), w'' = u(6)
% Here, INDEXSTARTDER = [1 4].

% First, we need to find out the maximum diffOrder of each variable appearing in
% the problem. We loop through the each component of the evaluation tree:
totalDiffOrders = zeros(1, numArgs);
for wCounter = 1:length(fevalResult)
    totalDiffOrders = max(totalDiffOrders, ...
        fevalResult(wCounter).tree.diffOrder);
end

% We always start indexing the first variable and its derivative(s) at 1. The
% index of the other variables depend on the cumulative diffOrders of the
% variables with lower indices.
indexStart = [ 1, cumsum(totalDiffOrders(1:end-1)) + 1 ];
% To get the indices of the derivatives as well, we shift all the previous
% indices right by 1, in cumulative fashion:
indexStartDer = indexStart + (0:numArgs-1);

% COEFFARG will be used to evaluate the functions that gives us information
% about the coefficients multiplying the highest order derivative in each
% equation. The vector has to be equally long to the total number of derivatives
% appearing in the problem; we'll then change one of the entries to 1 at a time
% to get the coefficient information.
coeffArg = zeros(1, indexStartDer(end) + totalDiffOrders(end));

% Go through each componenent from the result of evaluating FUNIN,
% and change it to infix format.
for wCounter = 1:length(fevalResult)
    
    % The current result we're looking at.
    res = fevalResult(wCounter);
    
    % Current diffOrders
    diffOrders = res.tree.diffOrder;
    
    % Ensure that we never have the highest derivatives of more than one
    % variable appearing in a single equation:
    if ( sum(totalDiffOrders == diffOrders) > 1 )
        error('CHEBFUN:TREEVAR:toFirstOrder:diffOrders', ...
            ['The highest order derivative of more than one variable ' ...
            'appears to be\npresent in the same equation. ' ...
            'Unable to convert to first order format.'])
    end
    
    % Expand the tree, so that PLUS rather than TIMES is sitting at the top of
    % it.
    expTree = treeVar.expandTree(res.tree, totalDiffOrders);
    
    % Split the tree into derivative part and non-derivative part.
    [newTree, derTree] = treeVar.splitTree(expTree, totalDiffOrders);
    
    % If newTree is empty, we only have a derivative part in the expression,
    % e.g. diff(u) = 0. We must replace it with a 0, as otherwise, we can't
    % evaluate the resulting odeFun in the ODE solvers.
    if ( isempty(newTree) )
        newTree = 0;
    end
    
    % Find what argument corresponds to the highest derivative one in the
    % current expression we're looking at. This will also be the order in which
    % we store the outputs from converting individual expressions to first order
    % form -- if the input is of the form @(x,u,v) [diff(v) - u; diff(u) - v],
    % the equations have to be sorted so that they'll be correctly converted.
    maxDerLoc = find(expTree.diffOrder == totalDiffOrders);
    
    % Convert the derivative part to infix form.
    % Indicate that we are converting a coeffFun
    isCoeffFun = true;
    [infixDer, varArrayDer] = ...
        treeVar.tree2infix(derTree, maxDerLoc, indexStartDer, isCoeffFun);
    
    % Convert the infix form of the expression that gives us the coefficient
    % multiplying the highest order derivative appearing in the expression to an
    % anonymous function we can evaluate:
    coeffFun = treeVar.toAnon(infixDer, varArrayDer);
    
    % Reset coeffArg for next evaluation:
    coeffArg = 0*coeffArg;
    
    % Replace one of the 0s in coeffFun with 1 so that we can evaluate COEFFFUN:
    if ( maxDerLoc == numArgs )
        % The last variable in the system currently appears in the highest order
        % derivate.
        coeffArg(end) = 1;
    else
        % The variable with index maxDerLoc+1 is the next variable we need to
        % start numbering at. So subtract 1 for the index of the highest
        % derivative we're currently interested in.
        coeffArg(indexStartDer(maxDerLoc+1) - 1) = 1;
    end
    
    % Evaluate the COEFFFUN to the coefficient!
    coeffs{maxDerLoc} = coeffFun(t, coeffArg);
    
    % Now work with the remaining syntax tree of the current expression of
    % interest. We need to negate the syntax tree as we're moving it to the
    % right-hand side. But if it already starts with a unary minus, we can
    % simply remove it rather than doing a double negation:
    % [TODO: Remove double UMINUS]
    newTree = struct('method', 'minus', 'numArgs', 2, ...
        'left', rhs{wCounter}, 'right', newTree);
    % Convert current expression to infix form:
    isCoeffFun = false;
    [infix, varArray] = ...
        treeVar.tree2infix(newTree, maxDerLoc, indexStart, isCoeffFun);
    % Store the infix form and the variables that appeared in the anonymous
    % function.
    systemInfix{maxDerLoc} = infix;
    varArrays{maxDerLoc} = varArray;
end

% Convert all the infix expressions, coefficients and variables stored to an
% anonymous function we can evaluate and use as the RHS of our ODE system:
funOut = treeVar.toRHS(systemInfix, varArrays, coeffs, ...
    indexStart, totalDiffOrders);

end
