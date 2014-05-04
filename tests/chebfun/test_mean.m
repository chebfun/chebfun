function pass = test_mean(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

%% Two arguments:
f = chebfun(@sin, pref);
g = chebfun(@cos, pref);
h = chebfun(@(x) .5*(sin(x) + cos(x)), pref);
pass(1) = normest(mean(f, g) - h) < 10*epslevel(h);

f = chebfun(@(x) [sin(x), cos(x)], pref);
g = chebfun(@(x) [cos(x), 1i*exp(x)], pref);
h = .5*(f+g);
pass(2) = normest(mean(f, g) - h) < 10*epslevel(h);

%% One argument:
f = chebfun(@sin, pref);
pass(3) = abs(mean(f)) < epslevel(f);
pass(4) = abs(mean(f.')) < epslevel(f);

f = chebfun(@(x) [sin(x), x], pref);
pass(5) = norm(mean(f), inf) < epslevel(f);
pass(6) = norm(mean(f.'), inf) < epslevel(f);

f = chebfun(@(x) [sin(x), x], [0, 6], pref);
pass(7) = norm(mean(f) - sum(f)/6, inf) < vscale(f)*epslevel(f);
pass(8) = norm(mean(f.') - sum(f).'/6, inf) < vscale(f)*epslevel(f);

%% singular function: a finite case

% define the domain:
dom = [-2 7];

op = @(x) sin(100*x)./((x-dom(1)).^0.5.*(x-dom(2)).^0.5);
f = chebfun(op, dom, 'exps', [-0.5 -0.5], 'splitting', 'on');
m = mean(f);
m_exact = -0.01273522016443600i;
pass(9) = abs(m - m_exact) < 1e1*get(f,'epslevel')*abs(m_exact);

%% singular function: an infinite case

% define the domain:
dom = [-2 7];

op = @(x) sin(100*x)./((x-dom(1)).^1.5.*(dom(2)-x).^0.5);
f = chebfun(op, dom, 'exps', [-1.5 -0.5], 'splitting', 'on');
m = mean(f);
pass(10) = ( isinf(m) );

%% singular function: a NaN case

% define the domain:
dom = [-2 7];

op = @(x) sin(98*x)./((x-dom(1)).^1.5.*(dom(2)-x).^1.5);
f = chebfun(op, dom, 'exps', [-1.5 -1.5], 'splitting', 'on');
m = mean(f);
pass(11) = ( isnan(m) );

%% Test for functions defined on unbounded domain:

% Functions on [-inf inf]:

% Set the domain:
dom = [-Inf 2 Inf];

op1 = @(x) x.^2.*exp(-x.^2);
op2 = @(x) (1-exp(-x.^2))./x.^2 + 2;
f = chebfun({op1 op2}, dom);
M = mean(f);
pass(12) = isnan(M);

% Function defined on [0 Inf]:

% Specify the domain: 
dom = [0 Inf];

op = @(x) 0.75+sin(10*x)./exp(x);
f = chebfun(op, dom, 'splitting', 'on');
M = mean(f);
pass(13) = abs(M - 0.75) < 10*epslevel(f).*vscale(f);

end

function out = normest(f, dom)

% Generate a few random points to use as test values.
seedRNG(6178);
if ( nargin == 1 )
    x = 2 * rand(100, 1) - 1;
else
    x = sum(dom) * rand(10, 1) - dom(1);
end

out = norm(feval(f, x), inf);

end
