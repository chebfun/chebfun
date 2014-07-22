function pass = test_basic1
%% Basic TREEVAR computations
u = treeVar();
v = cos(u);
w = sin(u);
t = v + w;
fprintf('\nExpression tree for t = cos(u) + sin(u):\n\n')
treeVar.printTree(t.tree)
treeVar.plotTree(t.tree)
pass = 1;
end