function pass = test_toFirstOrder()
%TEST_TOFIRSTORDER    Check that we get the correct output when changing to
%                     firt order format.
%% Setup
dom = [-1, 4];
x = chebfun(@(x) x, dom);
tol = 1e-14;
%% Simple test
% We have the expression 5*(diff(u, 2) + 3*u), so expect the first order system
% to be:
%   u'(1) = u(2)
%   u'(2) = -3*u(1)
problemNo = 1;
myfun = @(u) 5*(diff(u, 2) + 3*u);
rhs = 0;
[anonFun, idx, domOut, coeffs, diffOrders] = treeVar.toFirstOrder(myfun, rhs, dom);
pass(1, problemNo) = norm(anonFun(1,[2 1]) - [1;-6]) < tol;
pass(2, problemNo) = (idx == 1);
pass(3, problemNo) = all(domOut == dom);
pass(4, problemNo) = norm(coeffs{1} - 5) < tol;
pass(5, problemNo) = all( diffOrders == 2);
%% Try with variables in the anonymous function:
problemNo = 2;
alpha = 4;
myfun = @(u) 3.5*(diff(u, 2) + alpha*u);
rhs = 5;
% Expect the first order system to be:
%   u'(1) = u(2)
%   u'(2) = rhs/3.5 - alpha*u(1)
[anonFun, idx, domOut, coeffs, diffOrders] = treeVar.toFirstOrder(myfun, rhs, dom);
pass(1, problemNo) = norm(anonFun(1,[2 1]) - [1; 5/3.5-4*2]) < tol;
pass(2, problemNo) = (idx == 1);
pass(3, problemNo) = all(domOut == dom);
pass(4, problemNo) = norm(coeffs{1} - 3.5) < tol;
pass(5, problemNo) = all( diffOrders == 2);

%% Try with a CHEBFUN in the anonymous function:
problemNo = 3;
myfun = @(u) 5*((x+1).*diff(u, 2) + u);
rhs = -5;
[anonFun, idx, domOut, coeffs, diffOrders] = treeVar.toFirstOrder(myfun, rhs, dom);
% Expect this to be [u(2); (-5-5*u(1))/(5*(x(-.5)+1))] = [1; -4]
pass(1, problemNo) = norm( anonFun(-.5,[2 1]) - [1;-6] ) < tol;
pass(2, problemNo) = (idx == 1 );
pass(3, problemNo) = all(domOut == dom);
pass(4, problemNo) = norm(coeffs{1} - 5*(x+1)) < tol;
pass(5, problemNo) = all( diffOrders == 2);

%% Another expression with a CHEBFUN in the expression
problemNo = 4;
myfun = @(u) diff(u, 2) + 2*sin(x).*u;
rhs = cos(x);
[anonFun, idx, domOut, coeffs, diffOrders] = treeVar.toFirstOrder(myfun, rhs, dom);
% Expect this to be [u(2); cos(x(.5)) - 2*sin(x(.5))]
pass(1, problemNo) = norm( anonFun(.5,[2 1]) - [1; cos(x(.5))-2*(sin(x(.5))*2)]) < tol;
pass(2, problemNo) = (idx == 1);
pass(3, problemNo) = all(domOut == dom);
pass(4, problemNo) = ( norm(coeffs{1} - 1) < tol);
pass(5, problemNo) = all( diffOrders == 2);

%% Chebfun at the start
problemNo = 5;
myfun = @(u) 3*(x+2).*((x+1).*diff(u, 2) + u);
rhs = 0;
[anonFun, idx, domOut, coeffs, diffOrders] = treeVar.toFirstOrder(myfun, rhs, dom);
% Expect this to be [u(2); -u(1))/(x(-.5)+1)] = [1; -4]
pass(1, problemNo) = norm( anonFun(-.5,[2 2]) - [2; -4]) < tol;
pass(2, problemNo) = (idx == 1);
pass(3, problemNo) = all(domOut == dom);
pass(4, problemNo) = ( norm(coeffs{1} - 3*(x+2).*(x+1)) < tol);
pass(5, problemNo) = all( diffOrders == 2);
res =  anonFun(-.5, [2 2]);
pass(5) = all( res == [2;-4]);


%% Introduce breakpoints
problemNo = 6;
myfun = @(u) cos(x).*diff(u, 2) + 2*abs(sin(pi*x)).*u;
rhs = cos(2*x);
[anonFun, idx, domOut, coeffs, diffOrders] = treeVar.toFirstOrder(myfun, rhs, dom);
correctFun = @(x,u) [u(2); (cos(2*x)-2*abs(sin(pi*x)).*u(1))./cos(x)];
pass(1, problemNo) = norm(anonFun(.5,[2 1]) - correctFun(.5, [2 1])) < tol;
pass(2, problemNo) = (idx == 1);
pass(3, problemNo) = norm(domOut - (-1:4)) < tol;
pass(4, problemNo) = norm(coeffs{1} - cos(x)) < tol;
pass(5, problemNo) = all( diffOrders == 2);

%% Breakpoints, higher order derivatives
problemNo = 7;
myfun = @(u) cos(x).*diff(u, 4) + 2*abs(sin(pi*x)).*diff(u,2);
rhs = cos(2*x);
[anonFun, idx, domOut, coeffs, diffOrders] = treeVar.toFirstOrder(myfun, rhs, dom);
correctFun = @(x,u) [u(2); u(3); u(4); ...
    (cos(2*x)-2*abs(sin(pi*x)).*u(3))./cos(x)];
pass(1, problemNo) = norm(anonFun(.5,[2 1 3 4]) - correctFun(.5, [2 1 3 4])) < tol;
pass(2, problemNo) = (idx == 1);
pass(3, problemNo) = norm(domOut - (-1:4)) < tol;
pass(4, problemNo) = norm(coeffs{1} - cos(x)) < tol;
pass(5, problemNo) = all( diffOrders == 4);

%% Simple coupled system, first order
% We have the equations 
%   5*(diff(u) + 3*v) = tanh(x)
%   cos(x)*diff(v) + sin(x).*u = 2
% so expect the first order system to be:
%   u'(1) = tanh(x)/5 - 3*u(2)
%   u'(2) = (2 - sin(x)*u(1))/cos(x)
problemNo = 8;
myfun = @(x,u,v) [5*(diff(u) + 3*v); cos(x).*diff(v) + sin(x).*u];
rhs = [tanh(x); 2];
[anonFun, idx, domOut, coeffs, diffOrders] = treeVar.toFirstOrder(myfun, rhs, dom);
correctFun = @(x,u) [tanh(x)/5 - 3*u(2); (2-sin(x).*u(1))/cos(x)];
evalPt = [2.4 2.3];
pass(1, problemNo) = norm(anonFun(1, evalPt) - correctFun(1, evalPt)) < tol;
pass(2, problemNo) = all(idx == [1 2]);
pass(3, problemNo) = all(domOut == dom);
pass(4, problemNo) = norm(coeffs{1} - 5) + norm(coeffs{2} - cos(x)) < tol;
pass(5, problemNo) = all( diffOrders == [1 1]);


%% Simple coupled system, second order
% We have the equations 
%   5*(diff(u,2) + 3*v) = tanh(x)
%   cos(x)*diff(v,2) + sin(x).*u = 2
% so expect the first order system to be:
%   u'(1) = u(2)
%   u'(2) = tanh(x)/5 - 3*u(3)
%   u'(3) = u(4)
%   u'(4) = (2 - sin(x)*u(1))/cos(x)
problemNo = 9;
myfun = @(x,u,v) [5*(diff(u,2) + 3*v); cos(x).*diff(v,2) + sin(x).*u];
rhs = [tanh(x); 2];
[anonFun, idx, domOut, coeffs, diffOrders] = treeVar.toFirstOrder(myfun, rhs, dom);
correctFun = @(x,u) [u(2); tanh(x)/5 - 3*u(3); u(4); (2-sin(x).*u(1))/cos(x)];
evalPt = [2 1 2.4 2.3];
pass(1, problemNo) = norm(anonFun(1, evalPt) - correctFun(1, evalPt)) < tol;
pass(2, problemNo) = all(idx == [1 3]);
pass(3, problemNo) = all(domOut == dom);
pass(4, problemNo) = norm(coeffs{1} - 5) + norm(coeffs{2} - cos(x)) < tol;
pass(5, problemNo) = all( diffOrders == 2);

%% Coupled system, second order, breakpoints
% We have the equations 
%   5*(diff(u,2) + 3*v) = tanh(x)
%   cos(x)*diff(v,2) + abs(sin(pi*x)).*u = 2
% so expect the first order system to be:
%   u'(1) = u(2)
%   u'(2) = tanh(x)/5 - 3*u(3)
%   u'(3) = u(4)
%   u'(4) = (2 - sin(x)*u(1))/cos(x)
problemNo = 10;
myfun = @(x,u,v) [5*(diff(u,2) + 3*v); 
    cos(x).*diff(v,2) + abs(sin(pi*x)).*u];
rhs = [tanh(x); 2];
[anonFun, idx, domOut, coeffs, diffOrders] = treeVar.toFirstOrder(myfun, rhs, dom);
correctFun = @(x,u) [u(2); tanh(x)/5 - 3*u(3); u(4); (2-abs(sin(pi*x)).*u(1))/cos(x)];
evalPt = [2 1 2.4 2.3];
pass(1, problemNo) = norm(anonFun(1, evalPt) - correctFun(1, evalPt)) < tol;
pass(2, problemNo) = all(idx == [1 3]);
pass(3, problemNo) = norm(domOut - (-1:4)) < tol;
pass(4, problemNo) = norm(coeffs{1} - 5) + norm(coeffs{2} - cos(x)) < tol;
pass(5, problemNo) = all( diffOrders == 2);

%% Four variables, mixed derivatives, mixed order
% We have the equations
%   diff(y,3) + w + diff(u) = exp(x)
%   5*diff(w) + diff(v) = x.^2 
%   5*(diff(u,2) + 3*v) = tanh(x)
%   cos(x)*diff(v,2) + diff(u) + diff(y,2) = 2
%
% so expect the first order system to include the variables:
%   u(1) = u, u(2) = u', u(3) = v, u(4) = v'; 
% and the first order reformulation to be
%   u'(1) = u(2)
%   u'(2) = tanh(x)/5 - 3*u(3)
%   u'(3)
%   u'(2) = tanh(x)/5 - 3*u(3)
%   u'(3) = u(4)
%   u'(4) = (2 - sin(x)*u(1))/cos(x)
problemNo = 10;
myfun = @(x,u,v) [5*(diff(u,2) + 3*v); 
    cos(x).*diff(v,2) + abs(sin(pi*x)).*u];
rhs = [tanh(x); 2];
[anonFun, idx, domOut, coeffs, diffOrders] = treeVar.toFirstOrder(myfun, rhs, dom);
correctFun = @(x,u) [u(2); tanh(x)/5 - 3*u(3); u(4); (2-abs(sin(pi*x)).*u(1))/cos(x)];
evalPt = [2 1 2.4 2.3];
pass(1, problemNo) = norm(anonFun(1, evalPt) - correctFun(1, evalPt)) < tol;
pass(2, problemNo) = all(idx == [1 3]);
pass(3, problemNo) = norm(domOut - (-1:4)) < tol;
pass(4, problemNo) = norm(coeffs{1} - 5) + norm(coeffs{2} - cos(x)) < tol;
pass(5, problemNo) = all( diffOrders == 2);


%% Coupled systems -- Unsupported format, highest order derivatives in same eqn
myfun = @(x,u,v) [diff(u,2) + diff(v,2); diff(u) + sin(v)];
rhs = [1;2];
try
    treeVar.toFirstOrder(myfun, rhs, dom);
catch ME
    % The highest order derivatives of u and v appear in the same line -- this
    % should give us an error.
    errorPass(1) = strcmp(ME.identifier, 'CHEBFUN:TREEVAR:toFirstOrder:diffOrders');
end

%% Coupled systems -- Nonlinearity in highest order derivative
myfun = @(x,u,v) [diff(u,2).*diff(v); diff(v,2) + sin(u)];
rhs = [1;2];
try
    treeVar.toFirstOrder(myfun, rhs, dom);
catch ME
    % We're multiplying the highest order derivative by a variable, this should
    % give an error.
    errorPass(2) = strcmp(ME.identifier, 'CHEBFUN:TREEVAR:expandTree:nonlinearity');
end

%% Combine the information
pass = [all(pass) errorPass];
end
