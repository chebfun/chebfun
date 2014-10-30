function pass = test_expand
%% Expanding factored expressions
u = treeVar();
s3 = 5*(diff(u, 2) + 3*u);
fprintf('\n\nExpression tree for s3:\n\n')
treeVar.printTree(s3.tree);
figure
treeVar.plotTree(s3.tree)
title('Expression tree for s3 = 5*(diff(u,2) + 3*u))')
diffOrder = s3.tree.diffOrder;
fprintf('\n\nExpanded expression tree for s3:\n\n')
expTree = treeVar.expandTree(s3.tree, diffOrder);
treeVar.printTree(expTree);
figure
treeVar.plotTree(expTree)
title('Expanded expression tree for s3')
pass = 1;
end