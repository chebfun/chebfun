% Test vertical scale invariance of FOURTECH construction.

function pass = test_scaleInvariance(pref)

% Get preferences.
if ( nargin == 0 )
    pref = fourtech.techPref();
end

testclass = fourtech();

%% 
% Test CLASSICCHECK.
pref.happinessCheck = 'classic';

% Choose a test function and make a FOURTECH.
F = @(x) sin(10*pi*x);
f = testclass.make(F, [], pref);

% Require scale*FOURTECH(f) = FOURTECH(scale*f)).
scale = 2^300;
f1 = testclass.make(@(x) F(x)*scale, [], pref);
pass(1) = ~any(f.coeffs - f1.coeffs/scale);

% Require FOURTECH(f)/scale = FOURTECH(f/scale).
f2 = testclass.make(@(x) F(x)/scale, [], pref);
pass(2) = ~any(f.coeffs - f2.coeffs*scale);

%% 
% Test STRICTCHECK. (No STRICTCHECK for the moment for FOURTECH objects.)

end
