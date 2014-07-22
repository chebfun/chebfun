function pass = test_splitting
%% Splitting up the tree
u = treeVar();
s3 = 5*(diff(u, 2) + 3*u);
diffOrder = s3.tree.diffOrder;
expTree = treeVar.expandTree(s3.tree, diffOrder);

[newTree, derTree] = treeVar.splitTree(expTree, diffOrder);
fprintf('\n\nDerivative part of the expression tree for s3:\n\n')
treeVar.printTree(derTree);
treeVar.plotTree(derTree)
title('Derivative part of 5*(diff(u, 2) + 3*u)');
fprintf('\n\nRemaining part of the expression tree for s3:\n\n')
treeVar.printTree(newTree);
figure
treeVar.plotTree(newTree)
title('Remaining part of 5*(diff(u, 2) + 3*u)')
pass = 1;
end