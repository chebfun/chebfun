% Test vertical scale invariance of TRIGTECH construction.

function pass = test_scaleInvariance(pref)

% Get preferences.
if ( nargin == 0 )
    pref = trigtech.techPref();
end

testclass = trigtech();

%% 
% Test CLASSICCHECK.
pref.happinessCheck = 'classic';

% Choose a test function and make a TRIGTECH.
F = @(x) sin(10*pi*x);
f = testclass.make(F, [], pref);

% Require scale*TRIGTECH(f) = TRIGTECH(scale*f)).
scale = 2^300;
f1 = testclass.make(@(x) F(x)*scale, [], pref);
pass(1) = ~any(f.coeffs - f1.coeffs/scale);

% Require TRIGTECH(f)/scale = TRIGTECH(f/scale).
f2 = testclass.make(@(x) F(x)/scale, [], pref);
pass(2) = ~any(f.coeffs - f2.coeffs*scale);

%% 
% Test STRICTCHECK. (No STRICTCHECK for the moment for TRIGTECH objects.)

end
