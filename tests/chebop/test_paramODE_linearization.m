%% Simplest problem
A = chebop(@(x,u,c) u, [-1 1], @(x,u,c) u(0)-c);
f = chebfun(@exp);
u = A \ [f; 0];

%% Slightly more complicated
op = @(x,u,p) diff(u,2) + u;
bc = @(x,u,p) [ u(0)
                u(1) - 2
                feval(diff(u), 0) - p ];
A = chebop(op, [0 1], bc);
sol = A\0

%% Pass an initial guess to the problem above
op = @(x,u,p) diff(u,2) + u;
bc = @(x,u,p) [ u(0)
                u(1) - 2
                feval(diff(lam),0) - p ];

x = chebfun('x', [0 1])
init = [0*x; 0];

A = chebop(op, [0 1], bc, init);
sol = A\0

%% Example 1 (fine)
L = chebop(@(t,x,p) diff(x)-1,[-1 1]);
L.lbc = @(x,p) x;
L.rbc = @(x,p) x-p;
X = L\0;
[x, p] = deal(X);
plot(x)
p
% Clearly, the solution is p=2 and x(t) is a straight line from (-1,0) to (1,2).
% So far, Chebfun works fine. But now we try the same problem but using .bc
% syntax:
%% Example 2 (not working)
L = chebop(@(t,x,p) diff(x)-1,[-1 1]);
L.bc = @(t,x,p) [x(-1);x(1)-p];
X = L\0;
[x, p] = deal(X);
plot(x)
p
% Chebfun fails with a complicated error message, which I do not understand.
% However, if I slightly %change the problem so that the parameter p appears in
% the operator (replace -1 by -p in the first line) and in turn removing the p
% from the bc (replace -p by -2 in the second line), then Chebfun again works
% fine and gives the same x(t) and p=1 as expected:
%% Example 3 (fine)
L = chebop(@(t,x,p) diff(x)-p,[-1 1]);
L.bc = @(t,x,p) [x(-1);x(1)-2];
X = L\0;
[x, p] = deal(X);
plot(x)
p

%% If we get all the way here without errors, we're happy
% (Will do more proper pass checking once code's closer to working)
pass = 1 