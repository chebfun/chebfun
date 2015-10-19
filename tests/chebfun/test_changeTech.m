function pass = test_changeTech(pref)
% Test CHEBFUN/CHANGETECH.

if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e-14;

%% CHEBTECH TO TRIGTECH.
f = chebfun(@(x) cos(x), [0 2*pi]);
g = changeTech(f, @trigtech);
pass(1) = isequal(get(g.funs{1}, 'tech'), @trigtech);
pass(2) = norm(f - g, inf) < tol;

%% TRIGTECH TO CHEBTECH.
f = chebfun(@(x) cos(x), [0 2*pi], 'trig');
tech = pref.tech;
g = changeTech(f, tech);
pass(3) = isequal(get(g.funs{1}, 'tech'), tech);
pass(4) = norm(f - g, inf) < tol;

end
