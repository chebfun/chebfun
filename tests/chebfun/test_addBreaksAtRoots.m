function pass = test_addBreaksAtRoots(pref)

if ( nargin == 0 )
    pref = chebpref();
end

%% 
% Test that impulses are exactly zero at the new roots.

% Scalar:
f = chebfun(@(x) sin(x)-.5, pref);
g = addBreaksAtRoots(f);
pass(1) = g.pointValues(2) == 0;

% Array-valued:
f = chebfun(@(x) [sin(x), sin(x)-.5], pref);
g = addBreaksAtRoots(f);
pass(2) = g.pointValues(2,1) == 0 && g.pointValues(3,2) == 0;

% TODO: Add more tests.

end
