function pass = test_iszero(pref)

if ( nargin == 0 )
    pref = chebtech.pref();
end

f = chebtech.constructor(0, pref);
f.values = [0 1 0 ; 0 0 NaN];

pass(1) = all(iszero(f) == [1 0 0]);

f.values = [0 NaN 1];
pass(2) = all(iszero(f) == [1 0 0]);

f.values = [0 NaN 1]';
pass(3) = all(iszero(f) == 0);

f.values = zeros(3, 1);
pass(4) = all(iszero(f) == 1);

f.values = NaN;
pass(5) = all(iszero(f) == 0);

end