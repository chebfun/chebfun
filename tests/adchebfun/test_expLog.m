% Test file for ADCHEBFUN exp, log related functions

function pass = test_expLog

% List of trigonometric functions to test.
explogFunctions = {@exp, @expm1, @log, @log10, @log1p, @log2};

% Tolerance for Taylor testing
tol = 1e-2;

% Initialise vector with pass information
pass = zeros(3, numel(explogFunctions));

% Do the tests.
for k = 1:numel(explogFunctions)   
    
    % Call the valueTesting method, which also returns linearity information
    [err, lin] = adchebfun.valueTesting(explogFunctions{k});
    
    % First, check that the computed function values match what we expect
    pass(1, k) = ( err == 0 );
    
    % Call the taylorTesting method
    [order1, order2] = adchebfun.taylorTesting(explogFunctions{k});
    % We expect all elements of ORDER1 to be close to 1, and of ORDER2 to be
    % close to 2.
    pass(2, k) = ( max(abs(order1 - 1)) < tol ) & ...
        ( max(abs(order2 - 2)) < tol );
    
    % Check that we received the correct linearity information
    pass(3, k) = ( lin == 0 );
end

end