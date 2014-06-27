tic
mu = 1;
vdpFun = @(y) diff(y, 2) - mu.*(1-y.^2).*diff(y) + y;
x = chebfun(@(x) x);
u = vdpFun(treeChebfun(x));
[newTree, derTree] = treeChebfun.splitTree(u.tree, u.tree.diffOrder);
newTree = struct('method', 'uminus', 'numArgs', 1, 'center', newTree);
[infix, varCounter, varArray] = treeChebfun.tree2infix(newTree);
anonFun = treeChebfun.toAnon(infix,varArray);
[t,y]=ode113(anonFun,[0 20],[2 0]);
plot(t,y(:,1));
toc

tic
[t,y]=ode113(@vdp1,[0 20],[2 0]);
toc