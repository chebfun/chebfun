% Test file for ADCHEBFUN Bessel functions

function pass = test_bessel

% List of airy functions to test, with different number of inputs
besselj1 = @(f) besselj(1, f);
besselj25 = @(f) besselj(2.5, f);

besselFunctions = {besselj1, besselj25};

% Tolerance for Taylor testing
tol = 1e-2;

% Initialise vector with pass information
pass = zeros(2, numel(besselFunctions));

% Do the tests.
for k = 1:numel(besselFunctions)
    % First, check that the computed function values match what we expect
    pass(1, k) = ( adchebfun.valueTesting(besselFunctions{k}) == 0 );
    
    % Call the taylorTesting method
    [order1, order2] = adchebfun.taylorTesting(besselFunctions{k});
    % We expect all elements of ORDER1 to be close to 1, and of ORDER2 to be
    % close to 2.
    pass(2, k) = ( max(abs(order1 - 1)) < tol ) & ...
        ( max(abs(order2 - 2)) < tol );    
end

end