function pass = test_innerProduct(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

singPref = pref;
singPref.blowup = true;

% Seed for random number:
seedRNG(6178);

%% Functions on [-inf inf]:

% Set the domain:
dom = [-Inf Inf];

opf = @(x) 2-exp(-x.^2);
opg = @(x) exp(-x.^2);
f = unbndfun(opf, struct('domain', dom));
g = unbndfun(opg, struct('domain', dom));

I = innerProduct(f, g);
IExact = (sqrt(pi)*(4 - sqrt(2)))/2;
err = abs(I - IExact);
pass(1) = err < 2e7*max(get(f,'epslevel')*get(f,'vscale'), ...
    get(g,'epslevel')*get(g,'vscale'));

%% Functions on [a inf]:

% Set the domain:
dom = [1 Inf];

opf = @(x) x;
opg = @(x) exp(-x);
pref = chebfunpref();
f = unbndfun(opf, struct('domain', dom, 'exponents', [0 1]), pref);
g = unbndfun(opg, struct('domain', dom));
warning('off', 'CHEBFUN:UNBNDFUN:sum:slowDecay');
I = innerProduct(f, g);
warning('off', 'CHEBFUN:UNBNDFUN:sum:slowDecay');
IExact = 2*exp(-1);
err = abs(I - IExact);
pass(2) = err < 2e8*max(get(f,'epslevel')*get(f,'vscale'), ...
    get(g,'epslevel')*get(g,'vscale'));

%% Functions on [-inf b]:

% Set the domain:
dom = [-Inf -3*pi];

opf = @(x) 1./x;
opg = @(x) 2./x;
f = unbndfun(opf, struct('domain', dom));
g = unbndfun(opg, struct('domain', dom));
I = innerProduct(f, g);
IExact = 2/(3*pi);
err = abs(I - IExact);
pass(3) = err < 1e5*get(f,'epslevel')*get(f,'vscale');

end
