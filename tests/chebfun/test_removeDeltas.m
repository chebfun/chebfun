function pass = test_removeDeltas(pref)

% Get preferences:
if ( nargin < 1 )
    pref = chebfunpref();
end

% Test that the deltas are removed.
x = chebfun(@(x) x, [0 1], pref);
f = dirac(x - 0.2) + diff(dirac(x-1), 2);
thenorm = norm(removeDeltas(f), Inf); % Is infinite if deltas are present.
pass(1) = ( thenorm < Inf );

% Test that pointValues are changed.
f = dirac(x);
val = feval(removeDeltas(f), 0);
pass(2) = ( val == 0 );

% Test that removeDeltas is the identity for classicfuns.
f = sin(x);
g = removeDeltas(f);
pass(3) = isequal(f, g);

end
