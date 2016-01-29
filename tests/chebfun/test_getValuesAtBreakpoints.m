% Test file for @chebfun/test_getValuesAtBreakpoints.m.

function pass = test_getValuesAtBreakpoints(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

f = chebfun(@(x) x, [-1 0 .5 1]);
tol = eps;
vals = chebfun.getValuesAtBreakpoints(f.funs, f.domain);
err = abs(vals - f.domain.');
pass(1) = all(err < tol);
vals = chebfun.getValuesAtBreakpoints(f.funs, f.domain, {@(x) x, @(x) x, 1});
err = abs(vals - f.domain.');
pass(2) = all(err < tol);

pass(3) = all(chebfun.getValuesAtBreakpoints(f.funs, f.domain, ...
    @(x) x + 100*(x==0)) == [-1 100 .5 1].');

pass(4) = all(chebfun.getValuesAtBreakpoints(f.funs, f.domain, @sign) == ...
    [-1 0 1 1].');

%% piecewise smooth chebfun: smoothfun + singfun & splitting off:

% define the domain:
dom = [-2 -1 0 1];

op1 = @(x) sin(x);
op2 = @(x) 1./((1+x).^0.5);
op3 = @(x) x+1;
op = {op1, op2, op3};
f = chebfun(op, dom, 'exps', [0 0 -0.5 0 0 0]);
vals = chebfun.getValuesAtBreakpoints(f.funs);
vals_exact = [op1(dom(1)) Inf mean([op2(dom(3)) op3(dom(3))]) op3(dom(end))];
pass(5) = ( norm(vals([1, 3:4]) - vals_exact([1, 3:4]).', inf) ...
    < 1e1*eps*norm(vals_exact([1, 3:4]), inf) ) && ( vals(2) == Inf );

%% Test for function defined on unbounded domain:

% define the domain:
dom = [-2 0 Inf];

op1 = @(x) exp(x) - x;
op2 = @(x) 0.75+sin(10*x)./exp(x);
op = {op1, op2};

f = chebfun(op, dom);
vals = chebfun.getValuesAtBreakpoints(f.funs);
vals_exact = [op1(dom(1)); mean([op1(dom(2)) op2(dom(2))]); 0.75];

pass(6) = ( norm(vals - vals_exact, inf) < 1e3*eps*vscale(f) );


end
