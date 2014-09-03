function [funOut, indexStart, problemDom] = toFirstOrderSystem(funIn, rhs, domain)
% Independent variable on the domain
t = chebfun(@(t) t, domain);

% The first argument to funIn must be the independent time variable.
% So the number of treeVar arguments needed is one less:
numArgs = nargin(funIn) - 1;
args = cell(numArgs, 1);
argsVec = zeros(1, numArgs);

% Populate the args cell
for argCount = 1:numArgs
    argsVec(argCount) = 1;
    args{argCount} = treeVar(argsVec, domain);
    % Reset the index vector
    argsVec = 0*argsVec;
end

% Evaluate FUNIN with the TREEVAR arguments:
fevalResult = funIn(t, args{:});

% Initialize cells to store the infix forms of the expressions, the
% coefficients multiplying the highest order derivatives and any
% variables that appear in the anonymous functions:
systemInfix = cell(length(fevalResult),1);
coeffs = systemInfix;
varArrays = systemInfix;

% We want to return all potential breakpoints caused by piecewise
% coefficients in the problem. So loop through fevalResult, and
% unionize the breakpoints:
problemDom = fevalResult(1).domain;
for resCounter = 2:length(fevalResult)
    problemDom = union(problemDom, fevalResult(resCounter).domain);
end

% First look at all diffOrders to ensure we start with the correct
% indices. INDEXSTART denotes at which index we should start
% indexing each variable from. E.g., in the coupled system
%   [v'' + w; v + w'']
% we will have v = u(1), v' = u(2), w = u(3), w' = u(4), so
% INDEXSTART = [1, 3], since v and its derivatives starts getting
% counted at 1, and w and its derivatives start getting counted at
% 3.
%
% The vector INDEXSTARTDER is similar, but here, we also assign the
% highest order derivative of each variable it's own index. This is
% so that later on, we can correctly evaluate the coefficient
% multiplying the highest order derivative in problem. Thus, in the
% example above, we'd have
%   v = u(1), v' = u(2), v'' = u(3), w = u(4), w' = u(5), w'' = u(6)
% Here, INDEXSTARTDER = [1 4].
indexStart = zeros(1, numArgs);
indexStartDer = indexStart;
totalDiffOrders = indexStart;
for wCounter = 1:length(fevalResult)
    % We use cumsum() to look at how many derivatives have appeared
    % already in the problem.
    newIndex =  [1 (cumsum(fevalResult(wCounter).tree.diffOrder(1:end-1)) + (1:(numArgs-1)))];
    newIndexDer =  ...
        [1 (cumsum(fevalResult(wCounter).tree.diffOrder(1:end-1) + 1) + (1:(numArgs-1)))];
    indexStart = max(indexStart, newIndex);
    indexStartDer = max(indexStartDer, newIndexDer);
    totalDiffOrders = max(totalDiffOrders, fevalResult(wCounter).tree.diffOrder);
end

% COEFFARG will be used to evaluate the functions that gives us
% information about the coefficients multiplying the highest order
% derivative in each equation. The vector has to be equally long to
% the total number of derivatives appearing in the problem; we'll
% then change one of the entries to 1 at a time to get the
% coefficient information.
coeffArg = zeros(1, indexStartDer(end) + totalDiffOrders(end));

% Go through each componenent from the result of evaluating FUNIN,
% and change it to infix format.
for wCounter = 1:length(fevalResult)
    
    % The current result we're looking at.
    res = fevalResult(wCounter);
    % Current diffOrders
    diffOrders = res.tree.diffOrder;
    
    % Expand the tree, so that PLUS rather than TIMES is sitting at
    % the top of it.
    expTree = treeVar.expandTree(res.tree, totalDiffOrders);
    
    % Split the tree into derivative part and non-derivative part.
    [newTree, derTree] = treeVar.splitTree(expTree, diffOrders);
    
    % If newTree is empty, we only have a derivative part in the
    % expression, e.g. diff(u) = 0. We must replace it with a 0, as
    % otherwise, we can't evaluate the resulting odeFun in the ODE
    % solvers.
    if ( isempty(newTree) )
        newTree = 0;
    end
    
    % Convert the derivative part to infix form.
    [infixDer, dummy, varArrayDer] = ...
        treeVar.tree2infix(derTree, wCounter, indexStartDer);
    
    % Find what argument corresponds to the highest derivative one
    % in the current expression we're looking at:
    maxDerLoc = find(expTree.diffOrder == max(diffOrders));
    % Convert the infix form of the expression that gives us the
    % coefficient multiplying the highest order derivative appearing
    % in the expression to an anonymous function we can evaluate:
    coeffFun = treeVar.toAnon(infixDer, varArrayDer);
    
    % Reset coeffArg for next evaluation:
    coeffArg = 0*coeffArg;
    
    % Replace one of the 0s in coeffFun with 1 so that we can
    % evaluate COEFFFUN:
    if ( maxDerLoc == numArgs )
        % The last variable in the system currently appears in the
        % highest order derivate.
        coeffArg(end) = 1;
    else
        % The variable with index maxDerLoc+1 is the next variable
        % we need to start numbering at. So subtract 1 for the index
        % of the highest derivate we're currently interested in.
        coeffArg(indexStartDer(maxDerLoc+1) - 1) = 1;
    end
    
    % Evaluate the COEFFFUN to the coefficient!
    coeffs{wCounter} = coeffFun(t, coeffArg);
    
    % Now work with the remaining syntax tree of the current
    % expression of interest. We need to negate the syntax tree as
    % we're moving it to the right-hand side. But if it already
    % starts with a unary minus, we can simply remove it rather than
    % doing a double negation:
    % [TODO: Remove double UMINUS]
    newTree = struct('method', 'minus', 'numArgs', 2, ...
        'left', rhs{wCounter}, 'right', newTree);
    % Convert current expression to infix form:
    [infix, varCounter, varArray] = ...
        treeVar.tree2infix(newTree, wCounter, indexStart);
    % Store the infix form and the variables that appeared in the
    % anonymous function.
    systemInfix{wCounter} = infix;
    varArrays{wCounter} = varArray;
end

% Convert all the infix expressions, coefficients and variables
% stored to an anonymous function we can evaluate and use as the RHS
% of our ODE system:
funOut = treeVar.toRHS(systemInfix, varArrays, coeffs, ...
    indexStart, totalDiffOrders);

end