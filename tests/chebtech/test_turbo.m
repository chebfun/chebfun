% Test file for chebtech "turbo" option.

function pass = test_turbo(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebtech.techPref();
end

prefTurbo = pref;
prefTurbo.useTurbo = true;

for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else
        testclass = chebtech2();
    end

    % Check that the right number of coefficients were calculated.
    f_plain = testclass.make(@exp, [], pref);
    f_turbo = testclass.make(@exp, [], prefTurbo);
    pass(n, 1) = length(f_turbo) == 2*length(f_plain);

    % Make sure real-ness is respected.
    pass(n, 2) = isreal(f_turbo);

    % Do a rough check on the coefficients' absolute accuracy.
    rho = exp(abs(log(eps)) / length(f_plain))^(2/3);
    k = (0:1:(length(f_turbo) - 1)).';
    c_turbo = f_turbo.coeffs;
    c_exact = 2*besseli(k, 1);
    c_exact(1) = c_exact(1)/2;
    err = abs(c_turbo - c_exact);
    tol = 1e2*eps*rho.^(-k);
    pass(n, 3) = all(err < tol);

    % Try an array-valued example.
    f_plain = testclass.make(@(x) [exp(x) 1./(x + 5)], [], pref);
    f_turbo = testclass.make(@(x) [exp(x) 1./(x + 5)], [], prefTurbo);
    pass(n, 4) = length(f_turbo) == 2*length(f_plain);

    pass(n, 5) = isreal(f_turbo);

    rho = exp(abs(log(eps)) / length(f_plain))^(2/3);
    k = (0:1:(length(f_turbo) - 1)).';
    c_turbo = f_turbo.coeffs;
    c_ex1 = 2*besseli(k, 1);
    c_ex1(1) = c_ex1(1)/2;
    c_ex2 = (1/sqrt(6))*((-1).^k)./((5 + sqrt(24)).^k);
    c_ex2(1) = 1/(2*sqrt(6));
    c_exact = [c_ex1 c_ex2];
    err = abs(c_turbo - c_exact);
    tol = 1e2*eps*rho.^(-k);
    pass(n, 6) = all(all(bsxfun(@lt, err, tol)));

    % Check behavior in the presence of fixedLength.
    prefTurboFixed = prefTurbo;
    prefTurboFixed.fixedLength = 75;

    f_plain = testclass.make(@exp, [], pref);
    f_turbo = testclass.make(@exp, [], prefTurboFixed);
    pass(n, 7) = length(f_turbo) == 75;

    rho = exp(abs(log(eps)) / length(f_plain))^(2/3);
    k = (0:1:(length(f_turbo) - 1)).';
    c_turbo = f_turbo.coeffs;
    c_exact = 2*besseli(k, 1);
    c_exact(1) = c_exact(1)/2;
    err = abs(c_turbo - c_exact);
    tol = 1e2*eps*rho.^(-k);
    pass(n, 8) = all(err < tol);

    % Check behavior in the presence of fixedLength with an array-valued input.
    f_plain = testclass.make(@(x) [exp(x) 1./(x + 5)], [], pref);
    f_turbo = testclass.make(@(x) [exp(x) 1./(x + 5)], [], prefTurboFixed);
    pass(n, 9) = length(f_turbo) == 75;

    rho = exp(abs(log(eps)) / length(f_plain))^(2/3);
    k = (0:1:(length(f_turbo) - 1)).';
    c_turbo = f_turbo.coeffs;
    c_ex1 = 2*besseli(k, 1);
    c_ex1(1) = c_ex1(1)/2;
    c_ex2 = (1/sqrt(6))*((-1).^k)./((5 + sqrt(24)).^k);
    c_ex2(1) = 1/(2*sqrt(6));
    c_exact = [c_ex1 c_ex2];
    err = abs(c_turbo - c_exact);
    tol = 1e2*eps*rho.^(-k);
    pass(n, 10) = all(all(bsxfun(@lt, err, tol)));
end

end
