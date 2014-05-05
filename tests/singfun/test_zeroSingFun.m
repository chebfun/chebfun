% Test file for singfun/zeroSingFun.m

function pass = test_zeroSingFun(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

%%
% Construct the zero SINGFUN:
f = singfun.zeroSingFun();

% No exponents:
pass(1) = ~any(f.exponents);

% Trivial smooth part:
pass(2) = iszero(f.smoothPart);

end
