function pass = test_nonlinSysDamping_C1(pref)
% Test 2x2 system (sin/cos) where damping is required
%
% Asgeir Birkisson, May 2014.

if ( nargin == 0 )
    pref = cheboppref;
end

tol = 1e-10;

% Smooth domain:
d = [-pi pi];
A = chebop(@(x,u,v) [u - diff(v,2) + (1-x.^2).*u.^2; diff(u,2) + sin(v)],d);
A.lbc = @(u,v) [u-1; v+1];
A.rbc = @(u,v) [v-1/2; diff(v)];
x = chebfun('x',d);
f = [ 0*x ; 0*x ];

%% CHEBCOLLOC1
pref.discretization = @chebcolloc1;
[u1, u2, info] = solvebvp(A, f, pref);

% Want to check BCs as well.
bcFunLeft = chebfun(A.lbc(u1,u2));
bcFunRight = chebfun(A.rbc(u1,u2));

pass(1) = info.error < tol;
pass(2) = norm(bcFunLeft(d(1))) < tol && norm(bcFunRight(d(end))) < tol;

end


