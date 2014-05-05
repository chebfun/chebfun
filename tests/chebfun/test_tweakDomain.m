function pass = test_tweakDomain(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% Tweak left endpoint:
[df, dg, j, k] = test([-1+eps, 1], [-1, 1], pref);
pass(1) = all(df == dg) && j == 1 && k == 1;

% Tweak right endpoint:
[df, dg, j, k] = test([-1, 1], [-1, 1+eps], pref);
pass(2) = all(df == dg) && j == 2 && k == 2;

% Tweak a midpoint:
[df, dg, j, k] = test([-1, 0, 1], [-1, eps, 1], pref);
pass(3) = all(df == dg) && j == 2 && k == 2;

% Not all points are tweaked:
[df, dg, j, k] = test([-1, -.5, 0, .5, 1], [-1, eps, 1], pref);
pass(4) = all(dg == [-1, 0, 1]) && j == 3 && k == 2;

% Test not moving to an integer:
[df, dg, j, k] = test([-1, .5+eps, 1], [-1, .5-eps, 1], pref);
pass(5) = all(df == dg) && j == 2 && k == 2;

% Test multiple tweaks:
a = [-1-eps, -.5, -.2+eps, eps, .4, .7,  1];
b = [-1, -.2+eps, eps, .4+2*eps, .5, 1+eps];
[df, dg, j, k] = test(a, b, pref);
df2 = [-1, -.5, -.2+eps, eps, .4+eps, .7, 1];
dg2 = [-1, -.2+eps, eps, .4+eps, .5, 1];
pass(6) = all(df == df2) && all(dg == dg2) && ...
    all(j == [1 5 7]) && all(k == [1 4 6]);

end

function [df, dg, j, k] = test(df, dg, pref)
% Make dummy CHEBFUN objects with the given domains
f = chebfun(1, df, pref);
g = chebfun(1, dg, pref);
% Tweak their domains
[f, g, j, k] = tweakDomain(f, g);
df = f.domain;
dg = g.domain;

end
