% Test file for singfun/isempty.m

function pass = test_isequal(pref)

% Get preferences if not given
if ( nargin < 1 )
    pref = singfun.pref;
end

%%
% creat an empty SINGFUN
f = singfun;
g = singfun;
% following Matlab conventions, two empty objects are equal
pass(1) = isequal(f,g);

%%
% create a zero SINGFUN
f = singfun.zeroSingFun();
g = singfun.zeroSingFun();
pass(2) = isequal(f,g);

%%
% create a non-zero SINGFUN
f = singfun(@(x) 1./(1+x), [-1, 0] );
g = singfun(@(x) 1./(1+x), [-1, 0] );

% Test
pass(1) = isequal(f,g);
pass(2) = ~isempty(g);
pass(3) = ~isempty(h);
