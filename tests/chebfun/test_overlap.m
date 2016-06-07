% Test file for @chebfun/overlap.m.

function pass = test_overlap(pref)

% Grab some preferences:
if ( nargin == 0 )
    pref = chebfunpref();
end

% TODO: Test for array-valued CHEBFUNS and quasimatrices.

% Test empty input.
f = chebfun();
g = chebfun();
[f2, g2] = overlap(f, g);
pass(1) = isempty(f2) && isempty(g2);

% Check behavior for input chebfuns with differing domains.
f = chebfun(@sin, [-1 -.5 0 1], pref);
g = chebfun(@sin, [-2 0 0.5 2], pref);

try
    [f2, g2] = overlap(f, g);
    pass(2) = false;
catch ME
    pass(2) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:overlap:domains');
end

% Check behavior in the basic case.
f = chebfun(@sin, [-1 -.5 0 0.5 1], pref);
g = chebfun(@sin, [-1 0 0.5 1], pref);

[f2, g2] = overlap(f, g);
xx = linspace(-1, 1);
pass(3) = isequal(f2.domain, g2.domain) && ...
    norm(feval(f, xx) - feval(f2, xx), inf) < 10*eps*vscale(f) && ...
    norm(feval(g, xx) - feval(g2, xx), inf) < 10*eps*vscale(g);

% [TODO]: Are these the same tests as above?
% Check correct behavior for higher-order impulses.  [TODO]:  Use a function
% with real higher-order impulses instead of just faking them like this.
[f2, g2] = overlap(f, g);
xx = linspace(-1, 1);
pass(4) = isequal(f2.domain, g2.domain) && ...
    norm(feval(f, xx) - feval(f2, xx), inf) < 10*eps*vscale(f) && ...
    norm(feval(g, xx) - feval(g2, xx), inf) < 10*eps*vscale(g);

[g2, f2] = overlap(g, f);
xx = linspace(-1, 1);
pass(5) = isequal(f2.domain, g2.domain) && ...
    norm(feval(f, xx) - feval(f2, xx), inf) < 10*eps*vscale(f) && ...
    norm(feval(g, xx) - feval(g2, xx), inf) < 10*eps*vscale(g);

%% Test on singular function: piecewise smooth chebfun - splitting on.

% define the domain:
dom = [-2 7];
domCheck = [dom(1)+0.1 dom(2)-0.1];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

pow1 = -0.3;
pow2 = -0.5;
op1 = @(x) (x - dom(2)).^pow1.*sin(100*x);
op2 = @(x) (x - dom(2)).^pow2.*cos(300*x);
f = chebfun(op1, dom, 'exps', [0 pow1], 'splitting', 'on');
g = chebfun(op2, dom, 'exps', [0 pow2], 'splitting', 'on');
[fout, gout] = overlap(f,g);
vals_fout = feval(fout, x);
vals_gout = feval(gout, x);
vals_f = feval(op1, x);
vals_g = feval(op2, x);

check = zeros(1,4);
check(1) = all( fout.domain == gout.domain );
check(2) = all( fout.domain == unique([f.domain, g.domain]) );
check(3) = ( norm(vals_fout - vals_f, inf) < ...
    1e4*eps*norm(vals_fout, inf) );
check(4) = ( norm(vals_gout - vals_g, inf) < ...
    1e3*eps*norm(vals_gout, inf) );

pass(6) = all( check );


%% Test for function defined on unbounded domain:

% Function defined on [0 Inf]:

% Specify the domain: 
dom = [0 Inf];
domCheck = [0 100];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

% A decaying function:
opf = @(x) 0.75+sin(10*x)./exp(x);
f = chebfun(opf, dom, 'splitting', 'on');

% Blow-up function:
opg = @(x) x.*(5+exp(-x.^3))./(x - dom(1));
g = chebfun(opg, dom, 'exps', [-1 0]);

[fout, gout] = overlap(f, g);

vals_fout = feval(fout, x);
vals_gout = feval(gout, x);
vals_f = feval(opf, x);
vals_g = feval(opg, x);

check = zeros(1,4);
check(1) = all( fout.domain == gout.domain );
check(2) = all( fout.domain == unique([f.domain, g.domain]) );
check(3) = ( norm(vals_fout - vals_f, inf) < ...
    1e3*eps*norm(vals_fout, inf) );
check(4) = ( norm(vals_gout - vals_g, inf) < ...
    1e3*eps*norm(vals_gout, inf) );


pass(7) = all( check );

end
