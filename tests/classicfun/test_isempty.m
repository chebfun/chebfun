% Test file for @classicfun/isempty.m

function pass = test_isempty(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

singPref = pref;
singPref.blowup = true;

% Set a domain for BNDFUN.
data.domin = [-2 7];

%%
% Test an empty BNDFUN:
f = bndfun();
pass(1) = isempty(f);
    
%%
% Test an non-empty BNDFUN:
f = bndfun(@sin, data); 
pass(2) = ~isempty(f);

%%
% Test an non-empty array-valued BNDFUN:
f = bndfun(@(x) [sin(x), cos(x)], data);
pass(3) = ~isempty(f);
    
%%
% Test an array of BNDFUN objects:
f = [ bndfun(@sin, data), bndfun(@sin, data) ];
pass(4) = ~isempty(f);

%% Tests for UNBNDFUN:

% Functions on [-inf inf]:

% Set the domain:
data.domain = [-Inf Inf];

op = @(x) (1-exp(-x.^2))./x;
f = unbndfun(op, data);
pass(5) = ~isempty(f);

% Blow-up function:
op = @(x) x.^2.*(1-exp(-x.^2));
f = unbndfun(op, struct('domain', data.domain, 'exponents', [2 2]), singPref);
pass(6) = ~isempty(f);

%% Functions on [-inf b]:

% Set the domain:
data.domain = [-Inf -3*pi];

% Array-valued function:
op = @(x) [exp(x) x.*exp(x) (1-exp(x))./x];
f = unbndfun(op, data);
pass(7) = ~isempty(f);

end
