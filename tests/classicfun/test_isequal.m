% Test file for @classicfun/isequal.m

function pass = test_isequal(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

singPref = pref;
singPref.blowup = true;

% Set a domain for BNDFUN.
data.domain = [-2 7];

%% 
% Run a few basic tests for BNDFUN.
f = bndfun(@(x) sin(x), data, pref);
g = f;
pass(1) = isequal(f, g) && isequal(g, f);
    
g = bndfun(@(x) cos(x), data, pref);
pass(2) = ~isequal(f, g);
    
g = bndfun(@(x) [sin(x), cos(x)], data, pref);
pass(3) = ~isequal(f, g);
    
f = g;
pass(4) = isequal(f, g);
    
g = bndfun(@(x) [sin(x), exp(x)], data, pref);
pass(5) = ~isequal(f, g);
    
%% 
% Test on singular BNDFUN.
pow1 = -0.5;
pow2 = -0.6;
op1 = @(x) (x - data.domain(2)).^pow1.*sin(x);
op2 = @(x) (x - data.domain(2)).^pow2.*(cos(x).^2+1);
f = bndfun(op1, struct('domain', data.domain, 'exponents', [0 pow1]), pref);
g = bndfun(op2, struct('domain', data.domain, 'exponents', [0 pow2]), pref);
pass(6) = ~isequal(f, g);
    
%% Tests for UNBNDFUN:

% Functions on [-inf inf]:

% Set the domain:
data.domain = [-Inf Inf];

op = @(x) (1-exp(-x.^2))./x;
f = unbndfun(op, data);
pass(7) = isequal(f, f);

% Blow-up function:
op = @(x) x.^2.*(1-exp(-x.^2));
g = unbndfun(op, struct('domain', data.domain, 'exponents', [2 2]), singPref);
pass(8) = ~isequal(f, g);
pass(9) = ~isequal(g, f);

%% Functions on [-inf b]:

% Set the domain:
data.domain = [-Inf -3*pi];

% Array-valued function:
op = @(x) [exp(x) x.*exp(x) (1-exp(x))./x];
f = unbndfun(op, data);
pass(10) = isequal(f, f);

end
