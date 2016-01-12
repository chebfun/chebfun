function pass = test_followpath(~)
% Test the path following capabilities of Chebfun. This is not extensive, mainly
% checking nothing crashes in basic mode.

% Initialize a false pass vector
pass = zeros(1, 4);


%% Test Bratu problem
N = chebop(@(x,u,lam) diff(u,2) + lam*exp(u), [0 1]);
N.lbc = @(u,lam) u;
N.rbc = @(u,lam) u;
lam0 = 0.01;
% Call method, no plotting, no printing
try
    [u, lamvec] = followpath(N, lam0, 'maxstepno', 6); %#ok<*ASGLU>
    pass(1) = 1;
catch ME
    rethrow(ME)
end

% Call method, specifying more options, and more outputs
try
    [u, lamvec, mvec] = followpath(N, lam0, ...
        'measure', @(u) u(.5), 'maxstepno', 5);
    pass(2) = 1;
catch ME
    rethrow(ME)
end

%% Test Herceg problem, fix solution at left, vary slope at right
d = [0 1];
ep = 2^-2;
% Start by finding initial solution U on curve to be passed
N = chebop(@(x,u) -ep^2*diff(u,2)+(u.^2+u-.75).*(u.^2+u-3.75), d);
N.lbc = 0;
N.rbc = 0;
u0 = N\0;
% Trace solution curve:
N = chebop(@(x,u,lam) -ep^2*diff(u,2) + (u.^2+u-.75).*(u.^2+u-3.75), d);
lam0 = feval(diff(u0),d(2)); % Initial value for LAMBDA
N.lbc = @(u, lam) u;
N.rbc = @(u, lam) diff(u) - lam;
try
    [u, lamvec, mvec, lamfun, mfun] = followpath(N, lam0, 'maxstepno', 6, ...
        'uinit', u0, 'measure', @(u)u(1), 'stepmax', 1,  'stepmax', 0.1);
    pass(3) = 1;
catch ME
    rethrow(ME)
end

%% Herceg problem, continuation on perturbation parameter
ep = 2^-2;
H = chebop(@(x,u,lam) -lam*diff(u,2)+(u.^2+u-.75).*(u.^2+u-3.75), [0 1]);
lam0 = ep^2;
H.lbc = @(u, lam) u;
H.rbc = @(u, lam) u;
measure = @(u) norm(diff(u), inf);
try
    [u, lamvec, mvec] = followpath(H, lam0, ...
        'measure', measure, 'direction', -1, 'stopfun', @(u,lam) lam < 2e-2, ...
        'stepmax', .1);
    pass(4) = 1;
catch ME
    rethrow(ME)
end

end