% Test file for funcheb2/isequal.

function pass = test_isequal(pref)

% Get preferences.
if ( nargin < 1 )
    pref = funcheb.pref;
end

%%
% Run a few very straightforward tests.

f = funcheb2(@(x) sin(x), pref);
g = f;
pass(1) = isequal(f, g) && isequal(g, f);

g = funcheb2(@(x) cos(x), pref);
pass(2) = ~isequal(f, g);

g = funcheb2(@(x) [sin(x) cos(x)], pref);
pass(3) = ~isequal(f, g);

f = g;
pass(4) = isequal(f, g);

g = funcheb2(@(x) [sin(x) exp(x)], pref);
pass(5) = ~isequal(f, g);

end
