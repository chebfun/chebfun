% Test file for funcheb1/fliplr.

function pass = test_fliplr(pref)

% Get preferences.
if ( nargin < 1 )
    pref = funcheb.pref;
end

%%
% Conduct a few very straightforward tests.

f = funcheb1(@(x) sin(x), pref);
pass(1) = isequal(f, fliplr(f));

f = funcheb1(@(x) [sin(x) cos(x)], pref);
g = funcheb1(@(x) [cos(x) sin(x)], pref);
pass(2) = isequal(fliplr(f), g);

end
