% Test file for fun/isempty.m

function pass = test_isempty(varargin)

% Set a domain for BNDFUN.
dom = [-2 7];

%%
% Test an empty BNDFUN:
f = bndfun();
pass(1) = isempty(f);
    
%%
% Test an non-empty BNDFUN:
f = bndfun(@sin, dom); 
pass(2) = ~isempty(f);

%%
% Test an non-empty array-valued BNDFUN:
f = bndfun(@(x) [sin(x), cos(x)], dom);
pass(3) = ~isempty(f);
    
%%
% Test an array of BNDFUN objects:
f = [ bndfun(@sin, dom), bndfun(@sin, dom) ];
pass(4) = ~isempty(f);

%% Tests for UNBNDFUN:

% Functions on [-inf inf]:

% Set the domain:
dom = [-Inf Inf];

op = @(x) (1-exp(-x.^2))./x;
f = chebfun(op, dom);
pass(5) = ~isempty(f);

% Blow-up function:
op = @(x) x.^2.*(1-exp(-x.^2));
pref.singPrefs.exponents = [2 2];
f = chebfun(op, dom, pref); 
pass(6) = ~isempty(f);

%% Functions on [-inf b]:

% Set the domain:
dom = [-Inf -3*pi];

% Array-valued function:
op = @(x) [exp(x) x.*exp(x) (1-exp(x))./x];
f = chebfun(op, dom);
pass(7) = ~isempty(f);

end