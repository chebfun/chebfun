% Test vertical scale invariance of chebtech construction.

function pass = test_scaleInvariance(pref)

if ( nargin == 0 )
    pref = chebtech.techPref();
end

% Initialise pass matrix:
for n = 1:2
    
    % Test for both chebtech1 and chetech2
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end
    
    pref.happinessCheck = 'classic';
    
    % Choose a test function and make a chebtech:
    F = @(x) sin(10000*x);
    f = testclass.make(F, [], pref);
    
    % Require scale*chebtech(f) = chebtech(scale*f)):
    scale = 2^300;
    f1 = testclass.make(@(x) F(x)*scale, [], pref);
    pass(n,1) = ~any(f.coeffs - f1.coeffs/scale);

    % Require chebtech(f)/scale = chebtech(f/scale):
    f2 = testclass.make(@(x) F(x)/scale, [], pref);
    pass(n,2) = ~any(f.coeffs - f2.coeffs*scale);
    
    pref.happinessCheck = 'strict';
    
    % Choose a test function and make a chebtech:
    F = @(x) sin(10000*x);
    f = testclass.make(F, [], pref);

    % Require scale*chebtech(f) = chebtech(scale*f)):
    scale = 2^300;
    f1 = testclass.make(@(x) F(x)*scale, [], pref);
    pass(n,3) = ~any(f.coeffs - f1.coeffs/scale);

    % Require chebtech(f)/scale = chebtech(f/scale):
    f2 = testclass.make(@(x) F(x)/scale, [], pref);
    pass(n,4) = ~any(f.coeffs - f2.coeffs*scale);
    
end

end


