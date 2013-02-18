% Test file for funcheb2/fliplr.

function pass = test_fliplr(pref)

% Get preferences.
if ( nargin < 1 )
    pref = funcheb2.pref();
end

%%
% Conduct a few very straightforward tests.

f = funcheb2(@(x) sin(x), pref);
pass(1) = isequal(f, fliplr(f));

f = funcheb2(@(x) [sin(x) cos(x)], pref);
g = funcheb2(@(x) [cos(x) sin(x)], pref);
pass(2) = isequal(fliplr(f), g);

end
