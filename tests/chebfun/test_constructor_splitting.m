% Test file for chebfun constructor (splitting).

function pass = test_constructor_splitting(pref)

% Grab some preferences:
if ( nargin == 0 )
    pref = chebfunpref();
end

tol = 2000;

% Test SQRT(X) on [0 1]:
F1 = @sqrt;
f1 = chebfun(F1, [0, 1], pref, 'splitting', 'on', 'blowup', 'off');
xx1 = linspace(f1.domain(1)+eps, f1.domain(end)-eps, 100);
pass(1) = norm(feval(f1, xx1) - feval(F1, xx1), inf) < ...
    tol*max(f1.epslevel.*f1.vscale);

% Test SQRT(1-X) on [0 1]:
F2 = @(x) sqrt(1-x);
f2 = chebfun(F2, [0, 1], pref, 'splitting', 'on', 'blowup', 'off');
xx2 = linspace(f2.domain(1)+eps, f2.domain(end)-eps, 100);
pass(2) = norm(feval(f2, xx2) - feval(F2, xx2), inf) < ...
    tol*max(f2.epslevel.*f2.vscale);

% Test SQRT(1-X^2) on [-1 1]:
F3 = @(x) sqrt(1-x.^2);
f3 = chebfun(F3, [-1, 1], pref, 'splitting', 'on', 'blowup', 'off');
xx3 = linspace(f3.domain(1)+eps, f3.domain(end)-eps, 100);
pass(3) = norm(feval(f3, xx3) - feval(F3, xx3), inf) < ...
    tol*max(f3.epslevel.*f3.vscale);

% Test array-valued construction.  (Check for GitHub Issue #4.)
F4 = @(x) [sin(x) sign(x)];
f4 = chebfun(F4, [-1 1], pref, 'splitting', 'on', 'blowup', 'off');
xx4 = linspace(f4.domain(1)+eps, f4.domain(end)-eps, 100);
pass(4) = (norm(f4.domain - [-1 0 1], inf) < 10*eps) && ...
    (norm(feval(f4, xx4) - feval(F4, xx4), inf) < ...
    tol*max(f4.epslevel.*f4.vscale));

% Check for issue with call to merge on a function with multiple breakpoints.
F5 = @(x) sign(x - 0.1).*abs(x + 0.2).*sin(3*x);
f5 = chebfun(F5, [-1 1], pref, 'splitting', 'on', 'blowup', 'off');
xx5 = linspace(f5.domain(1)+eps, f5.domain(end)-eps, 100);
pass(5) = (norm(f5.domain - [-1 -0.2 0.1 1], inf) < 10*eps) && ...
    (norm(feval(f5, xx5) - feval(F5, xx5), inf) < ...
    tol*max(f5.epslevel.*f5.vscale));

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
pass(6) = norm(err, inf) < 1e1*epslevel(f)*vscale(f);

% % Test X*LOG(X) on [0 1]:
% F4 = @(x) x.*log(x);
% f4 = chebfun(F4, [0, 1], pref, 'splitting', 'on', 'blowup', 'off');
% xx4 = linspace(f4.domain(1)+eps, f4.domain(end)-eps, 100);
% pass(4) = norm(feval(f4, xx4) - feval(F4, xx4), inf) < ...
%     tol*max(f4.epslevel.*f4.vscale);
% 
% % Test (-X)*LOG(-X) on [-1 0]:
% F5 = @(x) (-x).*log(-x);
% f5 = chebfun(F5, [-1, 0], pref, 'splitting', 'on', 'blowup', 'off');
% xx5 = linspace(f5.domain(1)+eps, f5.domain(end)-eps, 100);
% pass(5) = norm(feval(f5, xx5) - feval(F5, xx5), inf) < ...
%     tol*max(f5.epslevel.*f5.vscale);
% 
% % Test (1-X)*LOG(1-X) on [0 1]:
% F6 = @(x) (1-x).*log(1-x);
% f6 = chebfun(F6, [0, 1], pref, 'splitting', 'on', 'blowup', 'off');
% xx6 = linspace(f6.domain(1)+eps, f6.domain(end)-eps, 100);
% pass(6) = norm(feval(f6, xx6) - feval(F6, xx6), inf) < ...
%     tol*max(f6.epslevel.*f6.vscale);
% pass(6) = pass(6) && numel(f6.funs) < 50;
% 
% % Test (1-X^2)*LOG(1-X^2) on [-1 1]:
% F7 = @(x) (1-x.^2).*log(1-x.^2);
% f7 = chebfun(F7, [-1, 1], pref, 'splitting', 'on', 'blowup', 'off');
% xx7 = linspace(f7.domain(1)+eps, f7.domain(end)-eps, 100);
% pass(7) = norm(feval(f7, xx7) - feval(F7, xx7), inf) < ...
%     tol*max(f7.epslevel.*f7.vscale);

end
