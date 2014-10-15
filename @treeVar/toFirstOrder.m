function [funOut, varIndex, problemDom] = toFirstOrder(funIn, rhs, domain)
% Independent variable on the domain
t = chebfun(@(t) t, domain);
arg = treeVar(1, domain);

% If funIn only has one input argument, we just give it a treeVar()
% argument. Otherwise, the first input will be the independent
% variable on the domain:
if ( nargin(funIn) == 1 )
    problemFun = funIn(arg);
else
    problemFun = funIn(t, arg);
end

% Return the potentially new breakpoints we obtained during the
% evaluation of problemFun:
problemDom = problemFun.domain;

% In the scalar case, varIndex will always be 1, as only one
% variable is involved:
varIndex = 1;

maxDifforder = problemFun.tree.diffOrder;

expTree = treeVar.expandTree(problemFun.tree, maxDifforder);

[newTree, derTree] = treeVar.splitTree(expTree, ...
    problemFun.tree.diffOrder);

% If newTree is empty, we only have a derivative part in the
% expression, e.g. diff(u) = 0. We must replace it with a 0, as
% otherwise, we can't evaluate the resulting odeFun in the ODE
% solvers.
if ( isempty(newTree) )
    newTree = 0;
end

[infixDer, dummy, varArrayDer] = treeVar.tree2infix(derTree, 1, 1);
coeffFun = treeVar.toAnon(infixDer, varArrayDer);
coeffArg = [zeros(1, expTree.diffOrder), 1];


coeff = {coeffFun(t, coeffArg)};

newTree = struct('method', 'minus', 'numArgs', 2, ...
    'left', rhs, 'right', newTree);
[infix, varCounter, varArray] = treeVar.tree2infix(newTree, 1, 1);
infix = {infix};
varArray = {varArray};
funOut = treeVar.toRHS(infix, varArray, coeff, 1, maxDifforder);
end
