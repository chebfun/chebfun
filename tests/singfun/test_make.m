% Test file for singfun/make.m

function pass = test_make(pref)

% Get preferences:
if ( nargin < 1 )
    pref = chebfunpref();
end

% tests on different calling sequences:

%%
% Providing only the function handle:
a = rand();
b = rand();
fh = @(x) sin(x)./((1+x).^a.*(1-x).^b);
f = singfun(fh);
g = f.make(fh);
pass(1) = isequal(f, g);

%%
% function handle and orders of singularities:
a = rand();
b = rand();
fh = @(x) sin(x).*(1+x).^a.*(1-x).^b;
data.exponents = [a, b];
f = singfun(fh, data);
data.exponents = [a, b];
g = f.make(fh, data);
pass(2) = isequal(f, g);

%%
% function handle and types of singularities:
a = ceil(10*rand);
b = ceil(10*rand);
fh = @(x) exp(x)./((1+x).^a.*(1-x).^b);
data.exponents = [];
data.singType = {'pole', 'pole'};
f = singfun(fh, data, pref);
data.exponents = [];
data.singType = {'pole', 'pole'};
g = f.make(fh, data, pref);
pass(3) = isequal(f, g);

%%
% function handle, orders of singularities, and preference:
a = rand();
b = rand();
fh = @(x) exp(sin(x))./((1+x).^a.*(1-x).^b);
data = struct();
data.exponents = [-a -b];
f = singfun(fh, data, pref);
data.exponents = [-a -b];
g = f.make(fh, data, pref);
pass(4) = isequal(f, g);

%%
% function handle, types of singularities, and preference: 
a = rand();
b = rand();
fh = @(x) sin(exp(cos(x))).*(1+x).^a.*(1-x).^b;
data.exponents = [];
data.singType = {'root', 'root'};
f = singfun(fh, data, pref);
data.exponents = [];
data.singType = {'root', 'root'};
g = f.make(fh, data, pref);
pass(5) = isequal(f, g);

%%
% all arguments are passed:
a = ceil(5*rand);
b = ceil(5*rand);
fh = @(x) exp(sin(x.^2))./((1+x).^a.*(1-x).^b);
data.exponents = [-a -b];
data.singType = {'pole', 'pole'};
f = singfun(fh, data, pref);
data.exponents = [-a -b];
data.singType = {'pole', 'pole'};
g = f.make(fh, data, pref);
pass(6) = isequal(f, g);
