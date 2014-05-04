% Test file for singfun/isfinite.m

function pass = test_isfinite(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

% The order of the exponents:
a = 0.64;
b = -0.64;
c = 1.28;
d = -1.28;

%%
% Check a few cases.

% fractional root at the left endpoint
f = singfun(@(x) (1 + x).^a.*exp(x), [a 0], {'root', 'none'}, [], [], pref);
pass(1) = isfinite(f);

% fractional pole at the left endpoint
f = singfun(@(x) (1 + x).^d.*sin(x), [d 0], {'sing', 'none'}, [], [], pref);
pass(2) = ~isfinite(f);

end
