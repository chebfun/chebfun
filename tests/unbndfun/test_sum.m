% Test file for unbndfun/sum.

function pass = test_sum(pref)

if ( nargin == 1 )
    pref = chebpref();
end

% Seed for random number:
seedRNG(6178);

%% Functions on [-inf inf]:

% Set the domain:
dom = [-Inf Inf];

op = @(x) exp(-x.^2);
f = unbndfun(op, dom);
I = sum(f);
IExact = sqrt(pi);
err = abs(I - IExact);
pass(1) = err < 5e2*get(f,'epslevel')*get(f,'vscale');

op = @(x) x.^2.*exp(-x.^2);
f = unbndfun(op, dom);
I = sum(f);
IExact = sqrt(pi)/2;
err = abs(I - IExact);
pass(2) = err < 2e5*get(f,'epslevel')*get(f,'vscale');

op = @(x) (1-exp(-x.^2))./x.^2;
f = unbndfun(op, dom);
I = sum(f);
IExact = 2*sqrt(pi);
err = abs(I - IExact);
pass(3) = err < 1e3*get(f,'epslevel')*get(f,'vscale');

%% Functions on [a inf]:

% Set the domain:
dom = [1 Inf];

op = @(x) exp(-x);
f = unbndfun(op, dom);
I = sum(f);
IExact = exp(-1);
err = abs(I - IExact);
pass(4) = err < 5e3*get(f,'epslevel')*get(f,'vscale');

op = @(x) x.*exp(-x);
f = unbndfun(op, dom);
I = sum(f);
IExact = 2*exp(-1);
err = abs(I - IExact);
pass(5) = err < 1e6*get(f,'epslevel')*get(f,'vscale');

op = @(x) (1-exp(-x))./x.^2;
f = unbndfun(op, dom);
I = sum(f);
IExact = 1 - exp(-1) - ei(-1);
err = abs(I - IExact);
pass(6) = err < 1e4*get(f,'epslevel')*get(f,'vscale');

op = @(x) 1./x.^2;
f = unbndfun(op, dom);
I = sum(f);
IExact = 1;
err = abs(I - IExact);
pass(7) = err < 1e4*get(f,'epslevel')*get(f,'vscale');

%% Functions on [-inf b]:

% Set the domain:
dom = [-Inf -3*pi];

op = @(x) exp(x);
f = unbndfun(op, dom);
I = sum(f);
IExact = exp(-3*pi);
err = abs(I - IExact);
pass(8) = err < 1e4*get(f,'epslevel')*get(f,'vscale');

op = @(x) x.*exp(x);
f = unbndfun(op, dom);
I = sum(f);
IExact = -exp(-3*pi)*(3*pi+1);
err = abs(I - IExact);
pass(9) = err < 1e4*get(f,'epslevel')*get(f,'vscale');

op = @(x) (1-exp(x))./x.^2;
f = unbndfun(op, dom);
I = sum(f);
IExact = (exp(-3*pi)*(exp(3*pi)-1))/(3*pi)-ei(-3*pi);
err = abs(I - IExact);
pass(10) = err < 1e5*get(f,'epslevel')*get(f,'vscale');

op = @(x) 1./x.^2;
f = unbndfun(op, dom);
I = sum(f);
IExact = 1/(3*pi);
err = abs(I - IExact);
pass(11) = err < 1e3*get(f,'epslevel')*get(f,'vscale');

end