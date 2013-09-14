function pass = test_addBreaksAtRoots(pref)

if ( nargin == 0 )
    pref = chebfun.pref();
end

%% 
% Test that impulses are exactly zero at the new roots.

% Scalar:
f = chebfun(@(x) sin(x)-.5, pref);
g = addBreaksAtRoots(f);
pass(1) = g.impulses(2) == 0;

% Array-valued:
f = chebfun(@(x) [sin(x), sin(x)-.5], pref);
g = addBreaksAtRoots(f);
pass(2) = g.impulses(2,1) == 0 && g.impulses(3,2) == 0;

% TODO: Add more tests.

end
