function pass = test_toFirstOrder()
%% Setup
dom = [-1, 4];
x = chebfun(@(x) x, dom);
%% Simple test
% We have the expression 5*(diff(u, 2) + 3*u), so expect the first order system
% to be:
%   u'(1) = u(2)
%   u'(2) = -3*u(1)
myfun = @(u) 5*(diff(u, 2) + 3*u);
rhs = 0;
anonFun = treeVar.toFirstOrder(myfun, rhs, dom);
pass(1) = all( anonFun(1,[2 1]) == [1 -6]');

%% Try with variables in the anonymous function:
alpha = 4;
myfun = @(u) 5*(diff(u, 2) + alpha*u);
rhs = 5;
anonFun = treeVar.toFirstOrder(myfun, rhs, dom);
pass(2) = all( anonFun(1,[2 1]) == [1;-7]);

%% Try with a CHEBFUN in the anonymous function:
myfun = @(u) 5*((x+1).*diff(u, 2) + u);
rhs = -5;
anonFun = treeVar.toFirstOrder(myfun, rhs, dom);
% Expect this to be [u(2); (-5-5*u(1))/(5*(x(-.5)+1))] = [1; -4]
pass(3) = all( anonFun(-.5, [2 1]) == [1;-6]);

%% Another expression with a CHEBFUN in the expression
myfun = @(u) diff(u, 2) + 2*sin(x).*u;
rhs = cos(x);
anonFun = treeVar.toFirstOrder(myfun, rhs, dom);
% Expect this to be [u(2); cos(x(.5)) - 2*sin(x(.5))] = [1; 1]
res = anonFun(.5,[2 1]);
pass(4) = norm(res - [1; cos(x(.5))-2*(sin(x(.5))*2)]) < 1e-14;

%% Chebfun at the start
myfun = @(u) 3*(x+2).*((x+1).*diff(u, 2) + u);
rhs = 0;
anonFun = treeVar.toFirstOrder(myfun, rhs, dom);
% Expect this to be [u(2); -u(1))/(x(-.5)+1)] = [1; -4]
res =  anonFun(-.5, [2 2]);
pass(5) = all( res == [2;-4]);

%% Coupled systems -- Unsupported format, highest order derivatives in same eqn
myfun = @(x,u,v) [diff(u,2) + diff(v,2); diff(u) + sin(v)];
rhs = [1;2];
try
    treeVar.toFirstOrder(myfun, rhs, dom);
catch ME
    % The highest order derivatives of u and v appear in the same line -- this
    % should give us an error.
    pass(6) = strcmp(ME.identifier, 'CHEBFUN:TREEVAR:toFirstOrder:diffOrders');
end

%% Coupled systems -- Nonlinearity in highest order derivative
myfun = @(x,u,v) [diff(u,2).*diff(v); diff(v,2) + sin(u)];
rhs = [1;2];
try
    treeVar.toFirstOrder(myfun, rhs, dom);
catch ME
    % We're multiplying the highest order derivative by a variable, this should
    % give an error.
    pass(7) = strcmp(ME.identifier, 'CHEBFUN:TREEVAR:expandTree:nonlinearity');
end
end