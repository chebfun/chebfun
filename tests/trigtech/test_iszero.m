% Test file for trigtech/iszero.m

function pass = test_iszero(pref)

% Get preferences.
if ( nargin < 1 )
    pref = trigtech.techPref();
end

testclass = trigtech();

f = testclass.make(0, pref);
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
