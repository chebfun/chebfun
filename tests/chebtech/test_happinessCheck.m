% Test file for chebtech/happinessCheck.m

function pass = test_happinessCheck(pref)

% Get preferences:
if ( nargin < 1 )
    pref = chebtech.pref;
end

pass = zeros(2, 6); % Pre-allocate pass matrix.
for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end

    % Set the tolerance:
    tol = 10*pref.chebtech.eps;
    
    %%
    % Test on a scalar-valued function:
    x = testclass.chebpts(33);
    f = @(x) sin(x);
    g = testclass.make(f(x));
    [ishappy, epslevel, tail] = happinessCheck(g, f, pref);
    pass(n, 1) = tail == 14;
    pass(n, 2) = ishappy && epslevel < tol;
    
    %%
    % Test on a vector-valued function:
    f = @(x) [sin(x) cos(x) exp(x)];
    g = testclass.make(f(x));
    [ishappy, epslevel, tail] = happinessCheck(g, f, pref);
    pass(n, 3) = abs(tail - 15) < 2;
    pass(n, 4) = ishappy && epslevel < tol;
    
    %%
    k = 32;
    m = k/2;
    x = testclass.chebpts(k+1);
    f = @(x) cos((2*k+m)*acos(x));
    
    % This should be happy, as aliasing fools the happiness test:
    pref.chebtech.sampletest = 0;
    g = testclass.make(f(x));
    [ishappy, epslevel, tail] = happinessCheck(g, f, pref);
    if (n == 1)
        pass(n, 5) = ( ishappy && tail == 15);
    else
        pass(n, 5) = ( ishappy && tail == 17);
    end
    
    % This should be unhappy, as sampletest fixes things:
    pref.chebtech.sampletest = 1;
    g = testclass.make(f(x));
    [ishappy, epslevel, tail] = happinessCheck(g, f, pref);
    pass(n, 6) = ~ishappy && tail == 33;
end

end
