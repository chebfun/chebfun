function pass = test_coeffs()
u = treeVar();
s3 = 5*(diff(u, 2) + 3*u);
diffOrder = s3.tree.diffOrder;
expTree = treeVar.expandTree(s3.tree, diffOrder);

[newTree, derTree] = treeVar.splitTree(expTree, diffOrder);
[infixDer, varArrayDer] = treeVar.tree2infix(derTree, 1, 1);
coeffFun = treeVar.toAnon(infixDer, varArrayDer);
coeffArg = [zeros(1, expTree.diffOrder), 1];
t = chebfun(@(t) t);
coeff = coeffFun(t, coeffArg);
pass = ( coeff == 5 );
end