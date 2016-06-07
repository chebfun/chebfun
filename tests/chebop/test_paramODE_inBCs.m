function pass = test_paramODE_inBCs(pref)
% Test solving a parameter dependent ODE, where the parameter appears in the
% BCs.
%
% See #1209.

if ( nargin == 0 )
    pref = cheboppref();
end
tol = 10*pref.bvpTol;

%% Second order problem with parameter in BC.
dom = [0 1];
x = chebfun(@(x) x, dom);
op = @(x,u,p) diff(u,2) + u;
bc = @(x,u,p) [ u(0)
                u(1) - 2
                feval(diff(u),0) - p ];
A = chebop(op, [0 1], bc);
[u,p] = deal(A\0);
res1 = A(x, u, p);
err(1) = norm(res1) + norm(bc(x,u,p));
%% Problems from mailing list (also included in #1209)
%% Example 1
L = chebop(@(t,x,p) diff(x)-1,[-1 1]);
L.lbc = @(x,p) x;
L.rbc = @(x,p) x-p;
U = L\0;
[u, p] = deal(U);
% The solution is p=2 and x(t) is a straight line from (-1,0) to (1,2). 
err(2) = norm(u-chebfun([0; 2])) + norm(p-2);

%% Example 2: Same above, but imposing conditions via the .bc field
L = chebop(@(t,x,p) diff(x)-1,[-1 1]);
L.bc = @(t,x,p) [x(-1);x(1)-p];
[u, p] = deal(L\0);
err(3) = norm(u-chebfun([0; 2])) + norm(p-2);
%% Example 3: Parameter appears in differential equation, rather than BCs
L = chebop(@(t,x,p) diff(x)-p,[-1 1]);
L.bc = @(t,x,p) [x(-1);x(1)-2];
[u, p] = deal(L\0);
err(4) = norm(u-chebfun([0; 2])) + norm(p-1);

%% Happy?
pass = err < tol;

end
