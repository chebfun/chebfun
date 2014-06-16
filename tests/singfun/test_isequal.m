% Test file for singfun/isequal.m

function pass = test_isequal(pref)

% Get preferences if not given
if ( nargin < 1 )
    pref = chebfunpref();
end

%%
% declare pass vector
%%
% create an empty SINGFUN
f = singfun();
g = singfun();
% following Matlab conventions, two empty objects are equal
pass(1) = isequal(f,g);

%%
% create a zero SINGFUN
f = singfun.zeroSingFun();
g = singfun.zeroSingFun();
pass(2) = isequal(f,g);

%%
% create a non-zero SINGFUN
data.exponents = [-1 0];
f = singfun(@(x) 1./(1+x), data);
g = singfun(@(x) 1./(1+x), data);

% Test
pass(3) = isequal(f,g);

%%
% create two different non-zero SINGFUN
data.exponents = [-1, 0];
data.singType = {'pole', 'none'};
f = singfun(@(x) 1./(1+x), data, pref);
data.exponents = [-1.8, 0];
data.singType = {'sing', 'none'};
g = singfun(@(x) 1./(1+x), data, pref);

% Test
pass(4) = ~isequal(f,g);

%%
% create two different non-zero SINGFUN
data.exponents = [-1, 0];
data.singType = {'pole', 'none'};
f = singfun(@(x) cos(x)./(1+x), data, pref);
data.exponents = [-1, 0];
data.singType = {'pole', 'none'};
g = singfun(@(x) sin(x)./(1+x), data, pref);

% Test
pass(5) = ~isequal(f,g);
