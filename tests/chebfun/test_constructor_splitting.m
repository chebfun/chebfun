% Test file for chebfun constructor (splitting).

function pass = test_constructor_splitting(pref)

% Grab some preferences:
if ( nargin == 0 )
    pref = chebfunpref();
end

seedRNG(6178);
tol = 1e9;

% Test SQRT(X) on [0 1]:
F1 = @sqrt;
f1 = chebfun(F1, [0, 1], pref, 'splitting', 'on', 'blowup', 'off');
xx1 = linspace(f1.domain(1)+eps, f1.domain(end)-eps, 100);
pass(1) = norm(feval(f1, xx1) - feval(F1, xx1), inf) < ...
    tol*max(eps*vscale(f1));

% Test SQRT(1-X) on [0 1]:
F2 = @(x) sqrt(1-x);
f2 = chebfun(F2, [0, 1], pref, 'splitting', 'on', 'blowup', 'off');
xx2 = linspace(f2.domain(1)+eps, f2.domain(end)-eps, 100);
pass(2) = norm(feval(f2, xx2) - feval(F2, xx2), inf) < ...
    tol*max(eps*vscale(f2));

% Test SQRT(1-X^2) on [-1 1]:
F3 = @(x) sqrt(1-x.^2);
f3 = chebfun(F3, [-1, 1], pref, 'splitting', 'on', 'blowup', 'off');
xx3 = linspace(f3.domain(1)+eps, f3.domain(end)-eps, 100);
pass(3) = norm(feval(f3, xx3) - feval(F3, xx3), inf) < ...
    tol*max(eps*vscale(f3));

% Test array-valued construction.  (Check for GitHub Issue #4.)
F4 = @(x) [sin(x) sign(x)];
f4 = chebfun(F4, [-1 1], pref, 'splitting', 'on', 'blowup', 'off');
xx4 = linspace(f4.domain(1)+eps, f4.domain(end)-eps, 100);
pass(4) = (norm(f4.domain - [-1 0 1], inf) < 10*eps) && ...
    (norm(feval(f4, xx4) - feval(F4, xx4), inf) < ...
    tol*max(eps*vscale(f4)));

% Check for issue with call to merge on a function with multiple breakpoints.
F5 = @(x) sign(x - 0.1).*abs(x + 0.2).*sin(3*x);
f5 = chebfun(F5, [-1 1], pref, 'splitting', 'on', 'blowup', 'off');
xx5 = linspace(f5.domain(1)+eps, f5.domain(end)-eps, 100);
pass(5) = (norm(f5.domain - [-1 -0.2 0.1 1], inf) < 10*eps) && ...
    (norm(feval(f5, xx5) - feval(F5, xx5), inf) < ...
    tol*max(eps*vscale(f5)));

%% Test a logical function:
f = chebfun(@(x) x > 0, [-1 1], 'splitting', 'on', pref);
x = chebfun('x', [-1 1], pref);
h = heaviside(x);
pass(6) = norm(f - h) < 10*eps;
    
% Test use of breakpoint detection in conjunction with construction from a cell
% array of function handles. (See issue #1151 on GitHub.)
f = chebfun({@(x) abs(x - 0.25), 0}, [0 0.5 1], 'splitting', 'on');
xx1 = linspace(0 + eps, 0.5 - eps, 20).';
err1 = norm(feval(f, xx1) - abs(xx1 - 0.25), inf);
xx2 = linspace(0.5 + eps, 1 - eps, 20).';
err2 = norm(feval(f, xx2), inf);
tol = 10*vscale(f)*eps;
pass(7) = max(err1, err2) < tol;

%% Test for function defined on unbounded domain:

% Function defined on [0 Inf]:

% Specify the domain: 
dom = [0 Inf];
domCheck = [0 100];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

op = @(x) 0.75+sin(10*x)./exp(x);
f = chebfun(op, dom, 'splitting', 'on');
fVals = feval(f, x);
fExact = op(x);
err = fVals - fExact;
pass(8) = norm(err, inf) < 1e2*eps*vscale(f);


%% Test SPLITTING ON with BLOWUP == 1:
op = @(x)tan(x);
f = chebfun(op, [-4 4], 'splitting', 'on', 'blowup', 1);

% Specify the domain: 
dom = [-4 -pi/2 pi/2 4];

% Generate random points to use as test values:
x1 = diff(dom(1:2)) * rand(100, 1) + dom(1);
x2 = diff(dom(2:3)) * rand(100, 1) + dom(2);
x3 = diff(dom(3:4)) * rand(100, 1) + dom(3);

err1 = op(x1)-f(x1);
err2 = op(x2)-f(x2);
err3 = op(x3)-f(x3);

err = [err1; err2; err3];
pass(9) = ( norm(err, inf) < 1e5*eps*vscale(f) );

%% Test for splitting on + unbndfun:

op = @(x)(sin(100*x)./exp(x.^2)+1).*(x.^2);
dom = [-inf inf];
dom_test = [-200 200];
x = diff(dom_test) * rand(100, 1) + dom_test(1);
f = chebfun (op, dom, 'exps', [2 2], 'splitting', 'on');
vals = f(x);
exact = op(x);
pass(10) = ( norm(vals-exact, inf) < 1e5*eps*vscale(f) );


% % Test X*LOG(X) on [0 1]:
% F4 = @(x) x.*log(x);
% f4 = chebfun(F4, [0, 1], pref, 'splitting', 'on', 'blowup', 'off');
% xx4 = linspace(f4.domain(1)+eps, f4.domain(end)-eps, 100);
% pass(4) = norm(feval(f4, xx4) - feval(F4, xx4), inf) < ...
%     tol*max(eps*vscale(f4));
% 
% % Test (-X)*LOG(-X) on [-1 0]:
% F5 = @(x) (-x).*log(-x);
% f5 = chebfun(F5, [-1, 0], pref, 'splitting', 'on', 'blowup', 'off');
% xx5 = linspace(f5.domain(1)+eps, f5.domain(end)-eps, 100);
% pass(5) = norm(feval(f5, xx5) - feval(F5, xx5), inf) < ...
%     tol*max(eps*vscale(f5));
% 
% % Test (1-X)*LOG(1-X) on [0 1]:
% F6 = @(x) (1-x).*log(1-x);
% f6 = chebfun(F6, [0, 1], pref, 'splitting', 'on', 'blowup', 'off');
% xx6 = linspace(f6.domain(1)+eps, f6.domain(end)-eps, 100);
% pass(6) = norm(feval(f6, xx6) - feval(F6, xx6), inf) < ...
%     tol*max(eps*vscale(f6));
% pass(6) = pass(6) && numel(f6.funs) < 50;
% 
% % Test (1-X^2)*LOG(1-X^2) on [-1 1]:
% F7 = @(x) (1-x.^2).*log(1-x.^2);
% f7 = chebfun(F7, [-1, 1], pref, 'splitting', 'on', 'blowup', 'off');
% xx7 = linspace(f7.domain(1)+eps, f7.domain(end)-eps, 100);
% pass(7) = norm(feval(f7, xx7) - feval(F7, xx7), inf) < ...
%     tol*max(eps*vscale(f7));

end
