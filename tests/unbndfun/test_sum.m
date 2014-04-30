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
tol = 5e3*get(f,'epslevel')*get(f,'vscale');
pass(1) = err < tol;

op = @(x) x.^2.*exp(-x.^2);
f = unbndfun(op, dom);
I = sum(f);
IExact = sqrt(pi)/2;
err = abs(I - IExact);
tol = 1e6*get(f,'epslevel')*get(f,'vscale');
pass(2) = err < tol;

op = @(x) (1-exp(-x.^2))./x.^2;
f = unbndfun(op, dom);
I = sum(f);
IExact = 2*sqrt(pi);
err = abs(I - IExact);
tol = 1e4*get(f,'epslevel')*get(f,'vscale');
pass(3) = err < tol;

% Blow-up function:
op = @(x) x.^2.*(1-exp(-x.^2));
pref.singPrefs.exponents = [2 2];
f = unbndfun(op, dom, [], [], pref); 
I = sum(f);
pass(4) = isequal(I, Inf);

% Blow-up function:
op = @(x) x;
pref.singPrefs.exponents = [1 1];
f = unbndfun(op, dom, [], [], pref); 
I = sum(f);
pass(5) = isnan(I);

%% Functions on [a inf]:

% Set the domain:
dom = [1 Inf];

op = @(x) exp(-x);
f = unbndfun(op, dom);
I = sum(f);
IExact = exp(-1);
err = abs(I - IExact);
tol = 1e4*get(f,'epslevel')*get(f,'vscale');
pass(6) = err < tol;

op = @(x) x.*exp(-x);
f = unbndfun(op, dom);
I = sum(f);
IExact = 2*exp(-1);
err = abs(I - IExact);
tol = 1e7*get(f,'epslevel')*get(f,'vscale');
pass(7) = err < tol;

op = @(x) (1-exp(-x))./x.^2;
f = unbndfun(op, dom);
I = sum(f);
% This exact value is obtained using Matlab's symbolic toolbox:
IExact = 0.851504493224078; 
err = abs(I - IExact);
tol = 1e4*get(f,'epslevel')*get(f,'vscale');
pass(8) = err < 2*tol;

op = @(x) 1./x.^2;
f = unbndfun(op, dom);
I = sum(f);
IExact = 1;
err = abs(I - IExact);
tol = 1e5*get(f,'epslevel')*get(f,'vscale');
pass(9) = err < tol;

% Blow-up function:
op = @(x) x.*(5+exp(-x.^3));
pref.singPrefs.exponents = [0 1];
f = unbndfun(op, dom, [], [], pref); 
I = sum(f);
pass(10) = isequal(I, Inf);

%% Functions on [-inf b]:

% Set the domain:
dom = [-Inf -3*pi];

op = @(x) exp(x);
f = unbndfun(op, dom);
I = sum(f);
IExact = exp(-3*pi);
err = abs(I - IExact);
tol = 1e4*get(f,'epslevel')*get(f,'vscale');
pass(11) = err < tol;

op = @(x) x.*exp(x);
f = unbndfun(op, dom);
I = sum(f);
IExact = -exp(-3*pi)*(3*pi+1);
err = abs(I - IExact);
tol = 1e4*get(f,'epslevel')*get(f,'vscale');
pass(12) = err < tol;

op = @(x) (1-exp(x))./x.^2;
f = unbndfun(op, dom);
I = sum(f);
% This exact value is obtained using Matlab's symbolic toolbox:
IExact = 0.106102535711326;
err = abs(I - IExact);
tol = 1e5*get(f,'epslevel')*get(f,'vscale');
pass(13) = err < tol;

op = @(x) 1./x.^2;
f = unbndfun(op, dom);
I = sum(f);
IExact = 1/(3*pi);
err = abs(I - IExact);
tol = 1e4*get(f,'epslevel')*get(f,'vscale');
pass(14) = err < tol;

op = @(x) 1./x.^2;
pref.singPrefs.exponents = [-2 0];
f = unbndfun(op, dom, [], [], pref);
I = sum(f);
IExact = 1/(3*pi);
err = abs(I - IExact);
tol = 1e1*get(f,'epslevel')*get(f,'vscale');
pass(15) = err < tol;

op = @(x) 0*x + 2;
f = unbndfun(op, dom);
I = sum(f);
pass(16) = isequal(I, Inf);

end
