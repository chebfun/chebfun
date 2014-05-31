function pass = test_plateauCheck(pref)

% TODO: Add more extensive testing.

if ( nargin == 0 )
    pref = chebfunpref();
end

pref.techPrefs.happinessCheck = @classicCheck;
f1 = chebtech2(@(x) [sin(x) cos(x)], [], [], pref);
pref.techPrefs.happinessCheck = @plateauCheck;
f2 = chebtech2(@(x) [sin(x) cos(x)], [], [], pref);
pass(1) = normest(f1 - f2) < 10*max(f2.epslevel);

end
