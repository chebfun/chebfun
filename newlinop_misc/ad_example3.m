% A simple example BVP.
clc, clear all, clear classes

% The 'CHEBOP':
N.op = @(x, u, v) [diff(u, 1) + v.^2 ; diff(v, 1) - u];


%%
% Initialise an ADCHEBFUN:
dom = [-1 0 1];
x = chebfun(@(x) x, dom);
u = adchebfun(@(u) 0*u, dom);
v = adchebfun(@(v) 0*v, dom);

u = seed(u, 1, 2);
v = seed(v, 2, 2);


%%
% Evaluate the operators to get a linearisation:

disp('nonlinear:')
w = u.^2;
isLinear = w.isConstant

' '
disp('linear:')
w = u + v;
isLinear = w.isConstant

' '
disp('nonlinear:')
w = sin(u);
isLinear = w.isConstant

' '
disp('linear:')
w = diff(u) + v;
isLinear = w.isConstant

' '
disp('nonlinear:')
w = diff(u).^2;
isLinear = w.isConstant

' '
disp('nonlinear:')
w = [diff(u, 2) + v ; diff(v) + u.^2];
isLinear = all([w.isConstant])