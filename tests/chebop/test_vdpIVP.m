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
toc
plot(t,y(:,1));


tic
[t,y]=ode113(@vdp1,[0 20],[2 0]);
toc

%% With treeVars
tic
mu = 25;
vdpFun = @(y) diff(y, 2) - mu.*(1-y.^2).*diff(y) + y;
u = vdpFun(treeVar());
[newTree, derTree] = treeVar.splitTree(u.tree, u.tree.diffOrder);
newTree = struct('method', 'uminus', 'numArgs', 1, 'center', newTree);
[infix, varCounter, varArray] = treeVar.tree2infix(newTree);
anonFun = treeVar.toAnon(infix,varArray);
[t,y]=ode113(anonFun,[0 40],[2 0]);
toc
plot(t,y(:,1));


tic
[t,y]=ode113(@vdp1,[0 20],[2 0]);
toc
figure, plot(t,y(:,1));
%% With chebfun.ode113
tic
mu = 1;
vdpFun = @(y) diff(y, 2) - mu.*(1-y.^2).*diff(y) + y;
u = vdpFun(treeVar());
[newTree, derTree] = treeVar.splitTree(u.tree, u.tree.diffOrder);
newTree = struct('method', 'uminus', 'numArgs', 1, 'center', newTree);
[infix, varCounter, varArray] = treeVar.tree2infix(newTree);
anonFun = treeVar.toAnon(infix,varArray);
y = chebfun.ode113(anonFun,[0 20],[2; 0])
toc
plot(y(:,1));

%% Chebop style
mu = 25;
vdpFun = @(y) diff(y, 2) - mu.*(1-y.^2).*diff(y) + y;
dom = [0 20];
N = chebop(vdpFun, dom);
N.lbc = @(u) [u - 2; diff(u)];
cheboppref.setDefaults('display','iter')
cheboppref.setDefaults('plotting','on')
[t, y] = N\0
plot(t,y(:,1));