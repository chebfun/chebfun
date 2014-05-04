function pass = test_addBreaksAtRoots(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

%% 
% Test that pointValues are exactly zero at the new roots.

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

pass(3) = ( norm(err, inf) < 1e2*epslevel(f)*norm(vals_check, inf) ) && ...
    ( norm( [dom(1); r_exact; dom(2)] - g.domain.', inf) < ...
    epslevel(f)*norm(r_exact, inf) );


%% Tests for functions defined on unbounded domain:

% Functions on [-inf inf]:

% Set the domain:
dom = [-Inf Inf];
domCheck = [-1e2 1e2];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

op = @(x) (1-exp(-x.^2))./x;
f = chebfun(op, dom);
g = addBreaksAtRoots(f);
rExact = 0;

vals_g = feval(g, x);
vals_check = feval(op, x);
err = vals_g - vals_check;
pass(4) = ( norm(err, inf) < 1e2*epslevel(f)*vscale(f) ) && ...
    ( norm( rExact - g.domain(2:end-1).', inf) < 2*epslevel(f)*vscale(f) );

% Blow-up function:
op = @(x) x.^2.*(1-exp(-x.^2))-2;
pref.singPrefs.exponents = [2 2];
f = chebfun(op, dom, pref); 
g = addBreaksAtRoots(f);
rExact = [-1.4962104914103104707 ; 1.4962104914103104707];

vals_g = feval(g, x);
vals_check = feval(op, x);
err = vals_g - vals_check;
pass(5) = ( norm(err, inf) < 5e3*epslevel(f)*vscale(f) ) && ...
    ( norm( rExact - g.domain(2:end-1).', inf) < epslevel(f)*vscale(f) );

% Functions on [a inf]:
dom = [0 Inf];
domCheck = [0 100];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

op = @(x) 0.15+sin(10*x)./exp(x);
f = chebfun(op, dom);
g = addBreaksAtRoots(f);
rExact = [0.33529141416564289113; 
          0.60061694515161002799;
          0.98375750309184861332;
          1.2042605667187311146;
          1.6619482204330474390;
          1.7760894757659030239];
      
vals_g = feval(g, x);
vals_check = feval(op, x);
err = vals_g - vals_check;
pass(6) = ( norm(err, inf) < epslevel(f)*vscale(f) ) && ...
    ( norm( rExact - g.domain(2:end-1).', inf) < 2*epslevel(f)*vscale(f) );

% TODO: Add more tests.

end
