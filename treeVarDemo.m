%% Basic TREEVAR computations
u = treeVar();
v = cos(u);
w = sin(u);
t = v + w;
fprintf('\nExpression tree for t:\n\n')
treeVar.printTree(t.tree)

%% Plot the tree
treeVar.plotTree(t.tree)
%% Introducing differentiation
x = chebfun(@(x) x);
myfun = @(u) 2 + diff(u,2);
s = myfun(u);
fprintf('\n\nExpression tree for s:\n\n')
treeVar.printTree(s.tree)
treeVar.plotTree(s.tree)
%% Nested differentiation
s2 = diff(diff(u)) + diff(u) + u;
fprintf('\n\nExpression tree for s2:\n\n')
treeVar.printTree(s2.tree)
treeVar.plotTree(s2.tree)
%% Expanding factored expressions:
s3 = 5*(diff(u, 2) + 3*u);
fprintf('\n\nExpression tree for s3:\n\n')
treeVar.printTree(s3.tree);
treeVar.plotTree(s3.tree)
diffOrder = s3.tree.diffOrder;
fprintf('\n\nExpanded expression tree for s3:\n\n')
expTree = treeVar.expandTree(s3.tree, diffOrder);
treeVar.printTree(expTree);
figure
treeVar.plotTree(expTree)
%% Splitting up the tree
[newTree, derTree] = treeVar.splitTree(expTree, diffOrder);
fprintf('\n\nDerivative part of the expression tree for s3:\n\n')
treeVar.printTree(derTree);
treeVar.plotTree(derTree)
fprintf('\n\nRemaining part of the expression tree for s3:\n\n')
treeVar.printTree(newTree);
figure
treeVar.plotTree(newTree)
%% Convert derTree to infix format to obtain coefficient:
[infixDer, dummy, varArrayDer] = treeVar.tree2infix(derTree, 1, 1);
infixDer %#ok<*NOPTS>
coeffFun = treeVar.toAnon(infixDer, varArrayDer);
coeffArg = [zeros(1, expTree.diffOrder), 1];
coeff = coeffFun(coeffArg);
%% Finish the conversion to a first order system:
newTree = struct('method', 'uminus', 'numArgs', 1, 'center', newTree);
[infix, varCounter, varArray] = treeVar.tree2infix(newTree, 1, 1);
funOut = treeVar.toRHS(infix, varArray, coeff)

%% Try evaluating funOut
% Recall that s3 = 5*(diff(u, 2) + 3*u), so expect the first order system to be:
%   u'(1) = u(2)
%   u'(2) = -3*u(1)
funOut(1, [2 1])
%% The method toFirstOrder takes care of the steps above:
dom = [0, 2];
x = chebfun(@(x) x, dom);
myfun = @(u) 5*(diff(u, 2) + 3*u);
anonFun = treeVar.toFirstOrder(myfun, dom);
anonFun(1,[2 1])

%%
% Try with variables in the anonymous function:
alpha = 4;
myfun = @(u) 5*(diff(u, 2) + alpha*u);
anonFun = treeVar.toFirstOrder(myfun, dom);
anonFun(1,[2 1]) % Expect this to be [u(2); -alpha*u(1)] = [1; -8]
%%
% And with a CHEBFUN in the anonymous function:
myfun = @(u) 5*((x+1).*diff(u, 2) + u);
anonFun = treeVar.toFirstOrder(myfun, dom);
anonFun(-.5,[2 1]) % Expect this to be [u(2); -u(1)/(x(-.5)+1)] = [1; -4]
%%
% And with a CHEBFUN in the anonymous function:
myfun = @(u) diff(u, 2) + x.*u;
anonFun = treeVar.toFirstOrder(myfun);
anonFun(-.5,[2 1]) % Expect this to be [u(2); -u(1)*(x(-.5)+1)] = [1; 1]


%% van der Pol
mu = 1;
vdpFun = @(y) diff(y, 2) - mu.*(1-y.^2).*diff(y) + y;
anonFun = treeVar.toFirstOrder(vdpFun);
tic
[t,y]=ode113(anonFun,[0 40],[2 0]);
fprintf('Solution time, conversion from 2nd order to 1st order: %4.4fs.\n', toc);
subplot(1,2,1), plot(t,y(:,1)); title('Automatically converted')

tic
[t,y]=ode113(@vdp1,[0 40],[2 0]);
fprintf('Solution time, built-in demo: %4.4fs.\n', toc);
subplot(1,2,2), plot(t,y(:,1)); title('Built in demo')

%% Crank up mu
mu = 40;
vdpFun = @(y) diff(y, 2) - mu.*(1-y.^2).*diff(y) + y;
anonFun = treeVar.toFirstOrder(vdpFun);
tic
[t,y]=ode113(anonFun,[0 40],[2 0]);
fprintf('Solution time, conversion 2nd->1st, larger mu: %4.4fs.\n', toc);
figure
plot(t,y(:,1)); title('Automatically converted')
%% With CHEBOPs
mu = 10;
dom = [0, 100];
vdpFun = @(y) diff(y, 2) - mu.*(1-y.^2).*diff(y) + y;
N = chebop(vdpFun, dom);
N.lbc = @(u) [u-2; diff(u)];
tic
y = N\0
fprintf('Solution time with chebops: %4.4fs.\n', toc);
figure
plot(y)
% Show where breakpoints got introduced
hold on, plot(y.domain, y(y.domain),'k*')
shg
%% We're resolving each piece well
hold off
plotcoeffs(y)

%% Chebfun computations with the solution:
plot(y)
hold on
r = roots(y-1.2);
plot(r, y(r), 'ro')
e = roots(diff(y));
plot(e, y(e), 'g*')