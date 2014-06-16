% Test file for singfun/isnan.m

function pass = test_isnan(pref)

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
data.exponents = [a 0];
data.singType = {'root', 'none'};
f = singfun(@(x) (1 + x).^a.*exp(x), data, pref);
pass(1) = ~isnan(f);

% fractional pole at the left endpoint
data.exponents = [d 0];
data.singType = {'sing', 'none'};
f = singfun(@(x) (1 + x).^d.*sin(x), data, pref);
pass(2) = ~isnan(f);

% fractional pole at the left endpoint
data.exponents = [d 0];
data.singType = {'sing', 'none'};
f = singfun(@(x) (1 + x).^d.*sin(x), data, pref);
g = NaN*f;
pass(3) = isnan(g);

end
