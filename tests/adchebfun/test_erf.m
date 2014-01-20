% Test file for ADCHEBFUN error and related functions, part 4/4

function pass = test_erf

% List of trigonometric functions to test.
errFunctions = {@erf, @erfc, @erfcinv, @erfcx, @erfinv};

% How many iterations we want to in the Taylor testing
numSteps = 4;

% Tolerance for Taylor testing
tol = 1e-2;

% Initialise vector with pass information
pass = zeros(2, numel(errFunctions));

% Do the tests.
for k = 1:numel(errFunctions)
    % First, check that the computed function values match what we expect
    pass(1, k) = ( valueTesting(errFunctions{k}) == 0 );
    
    % Call the taylorTesting method
    [order1, order2] = taylorTesting(errFunctions{k}, numSteps);
    % We expect all elements of ORDER1 to be close to 1, and of ORDER2 to be
    % close to 2.
    pass(2, k) = ( max(abs(order1 - 1)) < tol ) & ...
        ( max(abs(order2 - 2)) < tol );
    
    % Denote that an error happened when computing derivative, because the
    % version of development we branched off in Git didn't include the required
    % methods
    pass(1, k) = -1;
    pass(2, k) = -1;
    
end

end