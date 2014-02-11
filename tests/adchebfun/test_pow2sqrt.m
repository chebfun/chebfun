% Test file for ADCHEBFUN POW2() and SQRT()

function pass = test_pow2sqrt

% List of trigonometric functions to test.
pow2sqrtFunctions= {@pow2, @sqrt};

% How many iterations we want to in the Taylor testing
numSteps = 4;

% Tolerance for Taylor testing
tol = 1e-2;

% Initialise vector with pass information
pass = zeros(2, numel(pow2sqrtFunctions));

% Do the tests.
for k = 1:numel(pow2sqrtFunctions)
    % First, check that the computed function values match what we expect
    pass(1, k) = ( adchebfun.valueTesting(pow2sqrtFunctions{k}) == 0 );
    
    % Call the taylorTesting method
    [order1, order2] = adchebfun.taylorTesting(pow2sqrtFunctions{k}, numSteps);
    % We expect all elements of ORDER1 to be close to 1, and of ORDER2 to be
    % close to 2.
    pass(2, k) = ( max(abs(order1 - 1)) < tol ) & ...
        ( max(abs(order2 - 2)) < tol );
end

end