function pass = test_paramODE_linearization
% TEST_PARAMODE_LINEARIZATION    Test that we can correctly linearize and solve
% problems where unknown parameters appear in the boundary conditions.
tol = 1e-10;
%% Simplest problem
dom = [-1 1];
x = chebfun(@(x) x, dom);
A = chebop(@(x,u,c) diff(u), dom);
A.bc = @(x,u,c) [u(-1); u(0)- c];
f = chebfun(@exp);
[u, c] = deal(A \ f);
err = norm(A(u,c) - f) + norm(A.bc(x,u,c));
pass(1,1) = err < tol;
pass(2,1) = isa(u, 'chebfun') && isnumeric(c);

%% Same as above, but impose LBC as well
A = chebop(@(x,u,c) diff(u), dom);
A.bc = @(x,u,c) u(0)- c;
A.lbc = @(u,c) u;
[u, c] = deal(A \ f);
err = norm(A(u,c) - f) + norm(A.bc(x,u,c)) + norm(feval(A.lbc(u,c),-1));
pass(1,2) = err < tol;
pass(2,2) = isa(u, 'chebfun') && isnumeric(c);
%% Second order problem
dom = [0 1];
x = chebfun(@(x) x, dom);
op = @(x,u,p) diff(u,2) + u + p;
bc = @(x,u,p) [ u(0)
                u(1) - 2
                feval(diff(u),0) - p ];
A = chebop(op, dom, bc);
[u,p] = deal(A\0);
err = norm(A(u,p)) + norm(A.bc(x,u,p));
pass(1,3) = err < tol;
pass(2,3) = isa(u, 'chebfun') && isnumeric(p);


%% Original, with LBC and RBC
op = @(x,u,p) diff(u,2) + u + p;
lbc = @(u,p) [u; diff(u) - p];
rbc = @(u,p) u-2;
A = chebop(op, [0 1]);
A.lbc = lbc;
A.rbc = rbc;
[u, p] = deal(A \ 0);
err = norm(A(u,p)) + norm(feval(A.lbc(u,p), 0)) + norm(feval(A.rbc(u,p), 1));
pass(1,4) = err < tol;
pass(2,4) = isa(u, 'chebfun') && isnumeric(p);

%% Pass an initial guess to the problem above
op = @(x,u,p) diff(u,2) + u;
bc = @(x,u,p) [ u(0)
                u(1) - 2
                feval(diff(u),0) - p ];

x = chebfun(@(x) x, [0 1]);
init = [0*x; 0];

A = chebop(op, [0 1], bc, init);
[u,p] = deal(A\0);
err = norm(A(u,p)) + norm(A.bc(x,u,p));
pass(1,5) = err < tol;
pass(2,5) = isa(u, 'chebfun') && isnumeric(p);


%% Example 1 from Github issue
L = chebop(@(x,u,p) diff(u)-1, [-1 1]);
L.lbc = @(u,p) u;
L.rbc = @(u,p) u-p;
[u, p] = deal(L\0);
err = norm(L(u,p)) + norm(feval(L.lbc(u,p), -1)) + norm(feval(L.rbc(u,p), 1));
pass(1,6) = err < tol;
pass(2,6) = isa(u, 'chebfun') && isnumeric(p);

%% Example 2 from Github issue
L = chebop(@(x,u,p) diff(u)-1,[-1 1]);
L.bc = @(x,u,p) [u(-1);u(1)-p];
[u, p] = deal(L\0);
x = chebfun(@(x) x);
err = norm(L(u,p)) + norm(L.bc(x,u,p));
pass(1,7) = err < tol;
pass(2,7) = isa(u, 'chebfun') && isnumeric(p);
%% Example 3 from Github issue
L = chebop(@(x,u,p) diff(u)-p,[-1 1]);
L.bc = @(x,u,p) [u(-1);u(1)-2];
[u, p] = deal(L\0);
err = norm(L(u,p)) + norm(L.bc(x,u,p));
pass(1,8) = err < tol;
pass(2,8) = isa(u, 'chebfun') && isnumeric(p);

end