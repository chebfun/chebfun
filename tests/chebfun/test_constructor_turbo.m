function pass = test_constructor_turbo(pref)

% NB:  This test checks basic construction syntax only.  Accuracy testing is
% done in the corresponding test in tests/chebtech.

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

% Check that the right number of coefficients were calculated.
f_plain = chebfun(@exp, pref);
f_turbo = chebfun(@exp, 'turbo');
pass(1) = length(f_turbo) == 2*length(f_plain);

% Try an array-valued example.
f_plain = chebfun(@(x) [exp(x) 1./(x + 5)], pref);
f_turbo = chebfun(@(x) [exp(x) 1./(x + 5)], 'turbo');
pass(2) = length(f_turbo) == 2*length(f_plain);

% Try specifying a number of coefficients.
f = chebfun(@exp, 75, 'turbo');
pass(3) = length(f) == 75;

% Try specifying a number of coefficients for an array-valued input.
f = chebfun(@(x) [exp(x) 1./(x + 5)], 75, 'turbo');
pass(4) = length(f) == 75;

% Check that things work for constructions with breakpoints.
f_plain = chebfun(@exp, [-1 0 1], pref);
f_turbo = chebfun(@exp, [-1 0 1], 'turbo');
pass(5) = length(f_turbo) == 2*length(f_plain);

f_plain = chebfun(@(x) exp(x).*sign(x), 'splitting', 'on', pref);
f_turbo = chebfun(@(x) exp(x).*sign(x), 'splitting', 'on', 'turbo');
pass(6) = length(f_turbo) == 2*length(f_plain);

% Check that things work for constructions with exponents.
f_plain = chebfun(@(x) sin(x).*sqrt(1 + x), 'exps', [0.5 0], pref);
f_turbo = chebfun(@(x) sin(x).*sqrt(1 + x), 'exps', [0.5 0], 'turbo');
pass(7) = length(f_turbo) == 2*length(f_plain);

end
