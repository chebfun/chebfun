% Test file for ADCHEBFUN trigonometric and related functions, part 2/4

function pass = test_trig2

% List of trigonometric functions to test.
trigFunctions = {@acsch, @asec, @asecd, @asech, @asin, @asind, @asinh, @atan};

% Tolerance for Taylor testing
tol = 1e-2;

% Initialise vector with pass information
pass = zeros(2, numel(trigFunctions));

% Do the tests.
for k = 1:numel(trigFunctions)
    
    % First, check that the computed function values match what we expect
    pass(1, k) = ( adchebfun.valueTesting(trigFunctions{k}) == 0 );
    
    % Call the taylorTesting method
    [order1, order2] = adchebfun.taylorTesting(trigFunctions{k});
    % We expect all elements of ORDER1 to be close to 1, and of ORDER2 to be
    % close to 2.
    pass(2, k) = ( max(abs(order1 - 1)) < tol ) & ...
        ( max(abs(order2 - 2)) < tol );
end

end