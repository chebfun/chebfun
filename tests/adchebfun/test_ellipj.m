% Test file for ADCHEBFUN ellipj method

function pass = test_ellipj

% List of airy functions to test, with different inputs
ellipj075 = @(f) ellipj(f, 0.75);
ellipj1 = @(f) ellipj(f, 1);

ellipFunctions = {ellipj075, ellipj1};

% How many iterations we want to in the Taylor testing
numSteps = 6;

% Tolerance for Taylor testing
tol = 1e-2;

% Initialise vector with pass information
pass = zeros(2, numel(ellipFunctions));

% Do the tests.
for k = 1:numel(ellipFunctions)
    % First, check that the computed function values match what we expect
    pass(1, k) = ( adchebfun.valueTesting(ellipFunctions{k}, 3) == 0 );
    
    % Call the taylorTesting method
    [order1, order2] = adchebfun.taylorTesting(ellipFunctions{k}, numSteps, 3);
    % We expect all elements of ORDER1 to be close to 1, and of ORDER2 to be
    % close to 2.
    pass(2, k) = ( max(abs(order1 - 1)) < tol ) & ...
        ( max(abs(order2 - 2)) < tol );    
end

end