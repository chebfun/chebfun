% Test file for fun/isequal.m

function pass = test_isequal(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebpref();
end

% Set a domain for BNDFUN.
dom = [-2 7];

%% 
% Run a few basic tests for BNDFUN.
f = bndfun(@(x) sin(x), dom, [], [], pref);
g = f;
pass(1) = isequal(f, g) && isequal(g, f);
    
g = bndfun(@(x) cos(x), dom, [], [], pref);
pass(2) = ~isequal(f, g);
    
g = bndfun(@(x) [sin(x), cos(x)], dom, [], [], pref);
pass(3) = ~isequal(f, g);
    
f = g;
pass(4) = isequal(f, g);
    
g = bndfun(@(x) [sin(x), exp(x)], dom, [], [], pref);
pass(5) = ~isequal(f, g);
    
%% 
% Test on singular BNDFUN.
pow1 = -0.5;
pow2 = -0.6;
op1 = @(x) (x - dom(2)).^pow1.*sin(x);
op2 = @(x) (x - dom(2)).^pow2.*(cos(x).^2+1);
pref.singPrefs.exponents = [0 pow1];
f = bndfun(op1, dom, [], [], pref);
pref.singPrefs.exponents = [0 pow2];
g = bndfun(op2, dom, [], [], pref);
pass(6) = ~isequal(f, g);
    
%% Tests for UNBNDFUN:

% Functions on [-inf inf]:

% Set the domain:
dom = [-Inf Inf];

op = @(x) (1-exp(-x.^2))./x;
f = chebfun(op, dom);
pass(7) = isequal(f, f);

% Blow-up function:
op = @(x) x.^2.*(1-exp(-x.^2));
pref.singPrefs.exponents = [2 2];
g = chebfun(op, dom, pref); 
pass(8) = ~isequal(f, g);

%% Functions on [-inf b]:

% Set the domain:
dom = [-Inf -3*pi];

% Array-valued function:
op = @(x) [exp(x) x.*exp(x) (1-exp(x))./x];
f = chebfun(op, dom);
pass(9) = isequal(f, f);

end