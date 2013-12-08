function pass = test_addBreaksAtRoots(pref)

if ( nargin == 0 )
    pref = chebpref();
end

%% 
% Test that impulses are exactly zero at the new roots.

% Scalar:
f = chebfun(@(x) sin(x)-.5, pref);
g = addBreaksAtRoots(f);
pass(1) = g.impulses(2) == 0;

% Array-valued:
f = chebfun(@(x) [sin(x), sin(x)-.5], pref);
g = addBreaksAtRoots(f);
pass(2) = g.impulses(2,1) == 0 && g.impulses(3,2) == 0;

%% piecewise smooth chebfun: smoothfun + singfun & splitting on.

% define the domain:
dom = [-1 1];
domCheck = [dom(1)+0.1 dom(2)-0.1];

pow1 = -0.5;
pow2 = -1.2;
op = @(x) cos(100*x).*((x-dom(1)).^pow1).*((x-dom(2)).^pow2);
f = chebfun(op, dom, 'exps', [pow1 pow2], 'splitting', 'on');
g = addBreaksAtRoots(f);

% check values:

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

vals_g = feval(g, x);
vals_check = feval(op, x);
err = vals_g - vals_check;
pass(3) = norm(err-mean(err), inf) < 1e2*get(f,'epslevel')*norm(vals_check, inf);

% TODO: Add more tests.

end