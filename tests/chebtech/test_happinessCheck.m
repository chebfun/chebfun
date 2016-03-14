% Test file for chebtech/happinessCheck.m

function pass = test_happinessCheck(pref)

% Get preferences:
if ( nargin < 1 )
    pref = chebtech.techPref();
end

inputPref = pref;

for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end

    pref = inputPref;

    % Set the tolerance:
    tol = 10*pref.chebfuneps;
    
    %%
    % Test on a scalar-valued function:
    x = testclass.chebpts(33);
    f = @(x) sin(x);
    g = testclass.make(f(x));
    values = g.coeffs2vals(g.coeffs); 
    [ishappy, tail] = happinessCheck(g, f, values, [], pref);
    pass(n, 1) = tail == 14;
    pass(n, 2) = ishappy;
    
    %%
    % Test on an array-valued function:
    f = @(x) [sin(x) cos(x) exp(x)];
    g = testclass.make(f(x));
    values = g.coeffs2vals(g.coeffs); 
    [ishappy, tail] = happinessCheck(g, f, values, [], pref);
    pass(n, 3) = abs(tail - 15) < 2;
    pass(n, 4) = ishappy;
    
    %%
    k = 32;
    m = k/2;
    x = testclass.chebpts(k+1);
    f = @(x) cos((2*k+m)*acos(x));
    
    % This should be happy, as aliasing fools the happiness test:
    pref.sampleTest = 0;
    g = testclass.make(f(x));
    values = g.coeffs2vals(g.coeffs); 
    [ishappy, tail] = happinessCheck(g, f, values, [], pref);
    if (n == 1)
        pass(n, 5) = ( ishappy && tail == 15);
    else
        pass(n, 5) = ( ishappy && tail == 17);
    end
    
    % This should be unhappy, as sampletest fixes things:
    pref.sampleTest = 1;
    g = testclass.make(f(x));
    values = g.coeffs2vals(g.coeffs); 
    [ishappy, tail] = happinessCheck(g, f, values, [], pref);
    pass(n, 6) = ~ishappy && tail == 33;

    % g1 has a few coefficients that are small but not enough to satisfy
    % strictCheck, while g2 has just enough.  g1 will still pass with
    % classicCheck.
    pref.chebfuneps = 2^(-52);
    pref.happinessCheck = 'strict';
    f = @(x) sin(10*(x - 0.1));

    g1 = testclass.make(f(testclass.chebpts(39)));
    values1 = g1.coeffs2vals(g1.coeffs);
    ishappy1 = happinessCheck(g1, f, values1, [], pref);

    g2 = testclass.make(f(testclass.chebpts(41)));
    values2 = g2.coeffs2vals(g2.coeffs);
    ishappy2 = happinessCheck(g2, f, values2, [], pref);

    pref.happinessCheck = 'classic';
    ishappy3 = happinessCheck(g1, f, values1, [], pref);

    pass(n, 7) = ~ishappy1 && ishappy2 && ishappy3;

    % Test strictCheck with an array-valued input:
    pref.chebfuneps = 2^(-52);
    pref.happinessCheck = 'strict';
    f = @(x) [sin(10*(x - 0.1)) exp(x)];

    g1 = testclass.make(f(testclass.chebpts(39)));
    values1 = g1.coeffs2vals(g1.coeffs);
    ishappy1 = happinessCheck(g1, f, values1, [], pref);

    g2 = testclass.make(f(testclass.chebpts(41)));
    values2 = g2.coeffs2vals(g2.coeffs);
    ishappy2 = happinessCheck(g2, f, values2, [], pref);

    pass(n, 8) = ~ishappy1 && ishappy2;
    
    % Test plateauCheck with an array-valued input:
    p = pref;
    p.happinessCheck = @classicCheck;
    f1 = testclass.make(@(x) [sin(x) cos(x)], [], p);
    p.happinessCheck = @plateauCheck;
    f2 = testclass.make(@(x) [sin(x) cos(x)], [], p);
    pass(n, 9) = normest(f1 - f2) < 10*eps;

end

end
