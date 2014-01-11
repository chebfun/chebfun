function pass = test_constructor(pref)

if ( nargin == 1 )
    pref = chebpref();
end

%% Functions on [-inf inf]:
x = linspace(-10, 10);

F = @(x) exp(-x.^2);
f = unbndfun(F, [-inf, inf]);
pass(1) = norm(feval(f, x) - F(x), inf) < 10*get(f,'epslevel')*get(f,'vscale');

F = @(x) x.^2.*exp(-x.^2);
f = unbndfun(F, [-inf, inf]);
pass(2) = norm(feval(f, x) - F(x), inf) < 10*get(f,'epslevel')*get(f,'vscale');

F = @(x) (1-exp(-x.^2))./x;
f = unbndfun(F, [-inf, inf]); 
pass(3) = norm(feval(f, x) - F(x), inf) < 10*get(f,'epslevel')*get(f,'vscale');

%% Functions on [a inf]:
a = 1;
x = linspace(a, a + 20);

F = @(x) exp(-x);
f = unbndfun(F, [a, inf]);
pass(4) = norm(feval(f, x) - F(x), inf) < 10*get(f,'epslevel')*get(f,'vscale');

F = @(x) x.*exp(-x);
f = unbndfun(F, [a, inf]);
pass(5) = norm(feval(f, x) - F(x), inf) < 10*get(f,'epslevel')*get(f,'vscale');

F = @(x) (1-exp(-x))./x;
f = unbndfun(F, [a, inf]);
pass(6) = norm(feval(f, x) - F(x), inf) < 10*get(f,'epslevel')*get(f,'vscale');

F = @(x) 1./x;
f = unbndfun(F, [a, inf]);
pass(7) = norm(feval(f, x) - F(x), inf) < 10*get(f,'epslevel')*get(f,'vscale');

%% Functions on [-inf b]:
b = -pi;
x = linspace(b-20, b);

F = @(x) exp(x);
f = unbndfun(F, [-inf, b]);
pass(8) = norm(feval(f, x) - F(x), inf) < 10*get(f,'epslevel')*get(f,'vscale');

F = @(x) x.*exp(x);
f = unbndfun(F, [-inf, b]);
pass(9) = norm(feval(f, x) - F(x), inf) < 10*get(f,'epslevel')*get(f,'vscale');

F = @(x) (1-exp(x))./x;
f = unbndfun(F, [-inf, b]);
pass(10) = norm(feval(f, x) - F(x), inf) < 10*get(f,'epslevel')*get(f,'vscale');

F = @(x) 1./x;
f = unbndfun(F, [-inf, b]);
pass(11) = norm(feval(f, x) - F(x), inf) < 10*get(f,'epslevel')*get(f,'vscale');

%% MISC:

try
    f = unbndfun(@(x) exp(-x.^2), [0 1]);
    pass(12) = fail;
catch ME
    pass(12) = strcmp(ME.identifier, 'CHEBFUN:UNBNDFUN:BoundedDomain');
end
    
end
