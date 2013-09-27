% Test file for singfun/isequal.m

function pass = test_isequal(pref)

% Get preferences if not given
if ( nargin < 1 )
    pref = singfun.pref;
end

%%
% declare pass vector
pass = zeros(1,5);
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
f = singfun(@(x) 1./(1+x), [-1, 0] );
g = singfun(@(x) 1./(1+x), [-1, 0] );

% Test
pass(3) = isequal(f,g);

%%
% create two different non-zero SINGFUN
f = singfun(@(x) 1./(1+x), [-1, 0], {'pole', 'none'} );
g = singfun(@(x) 1./(1+x), [-1.8, 0], {'sing', 'none'} );

% Test
pass(4) = ~isequal(f,g);

%%
% create two different non-zero SINGFUN
f = singfun(@(x) cos(x)./(1+x), [-1, 0], {'pole', 'none'} );
g = singfun(@(x) sin(x)./(1+x), [-1, 0], {'pole', 'none'} );

% Test
pass(5) = ~isequal(f,g);