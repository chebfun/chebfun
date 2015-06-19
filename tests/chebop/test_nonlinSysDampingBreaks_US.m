function pass = test_nonlinSysDampingBreaks_US(pref)
% Test 2x2 system (sin/cos) where damping is required
%
% Asgeir Birkisson, May 2014.

if ( nargin == 0 )
    pref = cheboppref;
end

tol = 1e-10;

%% Piecewise:
d = [-pi 0 pi];
x = chebfun('x',d);
f = [ 0*x ; 0*x ];

%% ULTRAS:
pref.discretization = @ultraS;

A = chebop(@(x,u,v) [u - diff(v,2); diff(u,2) + cos(v)], d);
A.lbc = @(u,v) [u-1/2; v+1/4];
A.rbc = @(u,v) [v-1/4; u+1/2];

[u1, u2, info] = solvebvp(A, f, pref);

% Want to check BCs as well.
bcFunLeft = chebfun(A.lbc(u1,u2));
bcFunRight = chebfun(A.rbc(u1,u2));

% And check that we're continuous over breakpoint
u1jump = jump(u1, 0);
u2jump = jump(u2, 0);

pass(1) = abs(info.error) < tol;
pass(2) = norm(bcFunLeft(d(1))) < tol && norm(bcFunRight(d(end))) < tol;
pass(3) = norm(u1jump) < tol && norm(u2jump) < tol;

end
