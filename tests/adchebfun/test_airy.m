% Test file for ADCHEBFUN Airy function

function pass = test_airy

% List of airy functions to test, with different number of inputs
airy2 = @(f) airy(2, f);
funcList = {@airy, airy2};

% Tolerance for Taylor testing
tol = 1e-2;

% Initialise vector with pass information
pass = zeros(2, numel(funcList));

% Do the tests.
for k = 1:numel(funcList)
    % Call the valueTesting method, which also returns linearity information
    [err, lin] = adchebfun.valueTesting(funcList{k});
    
    % First, check that the computed function values match what we expect
    pass(1, k) = ( err == 0 );
    
    % Call the taylorTesting method
    [order1, order2] = adchebfun.taylorTesting(funcList{k});
    % We expect all elements of ORDER1 to be close to 1, and of ORDER2 to be
    % close to 2.
    pass(2, k) = ( max(abs(order1 - 1)) < tol ) & ...
        ( max(abs(order2 - 2)) < tol );
    
    % Check that we received the correct linearity information
    pass(3, k) = ( lin == 0 );
end

end