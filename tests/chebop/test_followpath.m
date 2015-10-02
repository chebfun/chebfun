function pass = test_followpath(~)
% Test the path following capabilities of Chebfun. This is not extensive, mainly
% checking nothing crashes in basic mode.

N = chebop(@(x,u,lam) diff(u,2) + lam*exp(u), [0 1]);
N.lbc = @(u,lam) u;
N.rbc = @(u,lam) u;
lam0 = 0.01;
% Call method, no plotting, no printing
try
    [u, lamvec] = followpath(N, lam0, 'maxstepno', 6);
    pass = 1;
catch ME
    pass = 0;
end