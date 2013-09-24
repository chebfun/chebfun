% A simple example BVP.
ccc

% The 'CHEBOP':
N.op = @(x, u, v) [diff(u, 1) + v ; diff(v, 1) - u];
N.lbc = @(u, v) u - 1;
N.bc = @(x, u, v) feval(v, .5);

%%
% Initialise an ADCHEBFUN:
x = chebfun('x', [-1 1]);
u = adchebfun(@(u) 0*u, [-1 1]);
u.jacobian = [linop.eye, linop.zeros];
v = adchebfun(@(v) 0*v, [-1 1]);
v.jacobian = [linop.zeros, linop.eye];

numVars = 2;

%%
% Evaluate the operators to get a linearisation:
w = N.op(x, u, v);
A = w(1).jacobian;
Ares = w(1).func;
for k = 2:numVars
    A = vertcat(A, w(k).jacobian);
    Ares = vertcat(Ares, w(k).func);
end

w = N.lbc(u, v);
B1 = w.jacobian;
B1res = w.func;
B1 = linop.evalAt(-1, [-1 1])*B1;
B1res = feval(B1res, -1);

w = N.bc(x, u, v);
B2 = w.jacobian;
B2res = w.func;

%%
% Assemble to a CHEBMATRIX:
% L = chebmatrix( { A } );
L = A;
L.constraints(1).op = B1;
L.constraints(1).value = -B1res;
L.constraints(2).op = B2;
L.constraints(2).value = -B2res;

% Make a CHEBMATRIX RHS:
rhs = chebmatrix( { x, x } );

%%
% Solve:
uv = L\rhs;

%%
% Display:
u = uv.blocks{1};
v = uv.blocks{2};
plot(u, 'b'); hold on
plot(v, 'r'), shg, hold off

u(-1) - 1
v(0.5)


