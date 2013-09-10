% Test file for @chebfun/isequal.m.

function pass = test_isequal(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfun.pref;
end

pref.chebfun.splitting = 1;

% Check empty case.
f = chebfun();
pass(1) = isequal(f, f);

% Check self equality.
f = chebfun(@(x) sin(x).*abs(x - 0.1), [-1 1], pref);
pass(2) = isequal(f, f);

% Check inequality for chebfuns with different row/column orientation.
pass(3) = ~isequal(f, f.');

% Check inequality for chebfuns of different dimensions.
g = chebfun(@(x) sin(x), [-1 1], pref);
pass(4) = ~isequal(f, g);

g = chebfun(@(x) [sin(x).*abs(x - 0.1) cos(x)], [-1 1], pref);
pass(5) = ~isequal(f, g);

% Check inequality for chebfuns with different domains.
g = chebfun(@(x) sin(x).*abs(x - 0.1), [-1+eps, 1+eps], pref);
pass(6) = ~isequal(f, g);

% Check inequality for chebfuns built from different functions.
g = chebfun(@(x) sin(x).*abs(x - 0.1) + eps, [-1 1], pref);
pass(7) = ~isequal(f, g);

end

