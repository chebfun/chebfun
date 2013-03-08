% Test file for funcheb1/poly.

function pass = test_poly(pref)

% Get preferences.
if ( nargin < 1 )
    pref = funcheb1.pref();
end

%%
% Check a few simple examples.

f = funcheb1(@(x) 3*ones(size(x)), pref);
p = poly(f);
pass(1) = (norm(p - 3, 'inf') < 10*f.epslevel);

f = funcheb1(@(x) 6.4*x - 3i, pref);
p = poly(f);
pass(2) = (norm(p - [6.4 (-3i)], 'inf') < 10*f.epslevel);

f = funcheb1(@(x) 2i*x.^5 - 3.2*x.^4 + 2*x.^2 - (1.2 + 3i), pref);
p = poly(f);
pass(3) = (norm(p - [2i (-3.2) 0 2 0 -(1.2 + 3i)], 'inf') < 10*f.epslevel);

%%
% Verify operation for vectorized funcheb1 objects.

f = funcheb1(@(x) [3*ones(size(x)), (6.4*x - 3i), (4*x.^2 - 2i*x + 3.7)], pref);
p = poly(f);
p_exact = [0 0     3;
           0 6.4   (-3i);
	   4 (-2i) 3.7];
pass(4) = (norm(p(:) - p_exact(:), 'inf') < 10*f.epslevel);

end
