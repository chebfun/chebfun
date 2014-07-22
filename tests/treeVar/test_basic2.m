function pass = test_basic2

%% Introducing differentiation
u = treeVar();
myfun = @(u) 2 + diff(u,2);
s = myfun(u);
fprintf('\n\nExpression tree for s = 2 + diff(u,2):\n\n')
treeVar.printTree(s.tree)
treeVar.plotTree(s.tree)
title(sprintf('Expression tree for s = 2 + diff(u,2)'))
%% Nested differentiation
s2 = diff(diff(u)) + diff(u) + u;
fprintf('\n\nExpression tree for s2 = diff(diff(u)) + diff(u) + u:\n\n')
treeVar.printTree(s2.tree)
figure
treeVar.plotTree(s2.tree)
title(sprintf(('Expression tree for s2 = diff(diff(u)) + diff(u) + u.')))

pass = 1;
end