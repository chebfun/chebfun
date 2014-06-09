% Test file for chebtech/iszero.m

function pass = test_iszero(pref)

if ( nargin == 0 )
    pref = chebtech.techPref();
end

for n = 1:2
    if ( n == 1)
        testclass = chebtech1();
    else
        testclass = chebtech2();
    end

    f = testclass.make(0, [], pref);
    f.coeffs = [0 1 0 ; 0 0 NaN];

    % iszero([0;0]) = 1, iszero([1;0]) = 0, iszero([0;NaN]) = 0.
    pass(n, 1) = all(iszero(f) == [1 0 0]);

    % iszero(0) = 1, iszero(NaN) = 0, iszero(1) = 0.
    f.coeffs = [0 NaN 1];
    pass(n, 2) = all(iszero(f) == [1 0 0]);
    f.coeffs = [0 NaN 1]';
    pass(n, 3) = all(iszero(f) == 0);

    % iszero(0) = 1.
    f.coeffs = zeros(3, 1);
    pass(n, 4) = all(iszero(f) == 1);

    % iszero(NaN) = 0.
    f.coeffs = NaN;
    pass(n, 5) = all(iszero(f) == 0);
end

end
