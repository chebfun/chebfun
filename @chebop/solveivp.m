function ty = solveivp(N, varargin)

disp('Initial value problem detected.')
% 
% problemFun = N.op(treeVar());
% [newTree, derTree] = treeVar.splitTree(problemFun.tree, ...
%     problemFun.tree.diffOrder);
% newTree = struct('method', 'uminus', 'numArgs', 1, 'center', newTree);
% [infix, varCounter, varArray] = treeVar.tree2infix(newTree);
% anonFun = treeVar.toAnon(infix,varArray);

anonFun = treeVar.toFirstOrder(N.op);
[t, y]=ode113(anonFun,N.domain,[2 0]);
ty.t = t;
ty.y = y;
end