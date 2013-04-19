function pass = test_iszero(pref)

if ( nargin == 0 )
    pref = chebtech.pref();
end

f = chebtech.constructor(0, pref);
f.values = [0 1 0 ; 0 0 NaN];

% iszero([0;0]) = 1, iszero([1;0]) = 0, iszero([0;NaN]) = 0.
pass(1) = all(iszero(f) == [1 0 0]);

% iszero(0) = 1, iszero(NaN) = 0, iszero(1) = 0.
f.values = [0 NaN 1];
pass(2) = all(iszero(f) == [1 0 0]);
f.values = [0 NaN 1]';
pass(3) = all(iszero(f) == 0);

% iszero(0) = 1.
f.values = zeros(3, 1);
pass(4) = all(iszero(f) == 1);

% iszero(NaN) = 0.
f.values = NaN;
pass(5) = all(iszero(f) == 0);

end