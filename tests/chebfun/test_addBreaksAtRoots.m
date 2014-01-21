function pass = test_addBreaksAtRoots(pref)

if ( nargin == 0 )
    pref = chebpref();
end

%% 
% Test that impulses are exactly zero at the new roots.

% Scalar:
f = chebfun(@(x) sin(x)-.5, pref);
g = addBreaksAtRoots(f);
pass(1) = g.pointValues(2) == 0;

% Array-valued:
f = chebfun(@(x) [sin(x), sin(x)-.5], pref);
g = addBreaksAtRoots(f);
pass(2) = g.pointValues(2,1) == 0 && g.pointValues(3,2) == 0;

%% piecewise smooth chebfun: smoothfun + singfun & splitting on.

% define the domain:
dom = [-2 7];
domCheck = [dom(1)+0.1 dom(2)-0.1];

pow1 = -0.5;
pow2 = -1.2;
op = @(x) cos(30*x).*((x-dom(1)).^pow1).*((x-dom(2)).^pow2);
f = chebfun(op, dom, 'exps', [pow1 pow2], 'splitting', 'off');
g = addBreaksAtRoots(f);

% check values:

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

vals_g = feval(g, x);
vals_check = feval(op, x);
err = vals_g - vals_check;

r_exact = (((-19:66)+1/2)*pi/30).';

pass(3) = ( norm(err-mean(err), inf) < ...
    1e2*get(f,'epslevel')*norm(vals_check, inf) ) && ...
    ( norm( [dom(1); r_exact; dom(2)] - g.domain.', inf) < ...
    get(f,'epslevel')*norm(r_exact, inf) );

% TODO: Add more tests.

end
