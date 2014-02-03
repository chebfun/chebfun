% Test file for ADCHEBFUN airy function

function pass = test_airy

% List of airy functions to test, with different number of inputs
airy2 = @(f) airy(2, f);
airyFunctions = {@airy, airy2};

% How many iterations we want to in the Taylor testing
numSteps = 6;

% Tolerance for Taylor testing
tol = 1e-2;

% Initialise vector with pass information
pass = zeros(2, numel(airyFunctions));

% Do the tests.
for k = 1:numel(airyFunctions)
    % First, check that the computed function values match what we expect
    pass(1, k) = ( adchebfun.valueTesting(airyFunctions{k}) == 0 );
    
    % Call the taylorTesting method
    [order1, order2] = adchebfun.taylorTesting(airyFunctions{k}, numSteps);
    % We expect all elements of ORDER1 to be close to 1, and of ORDER2 to be
    % close to 2.
    pass(2, k) = ( max(abs(order1 - 1)) < tol ) & ...
        ( max(abs(order2 - 2)) < tol );    
end

end