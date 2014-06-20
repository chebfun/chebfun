% Test file for ADCHEBFUN FRED

function pass = test_fred

K = @(s, t)  exp(-(s-t).^2);

% List of trigonometric functions to test.
func = @(u) fred(K, u);

% Tolerance for Taylor testing
tolOrder = 1e-2;
tolDiff = 1e-12;
% Initialise vector with pass information
pass = zeros(1, numel(func));

% Call the valueTesting method, which also returns linearity information
[err, lin] = adchebfun.valueTesting(func);

% First, check that the computed function values match what we expect
pass(1) = ( err == 0 );

% Call the taylorTesting method
[order1, order2, nDiff2] = adchebfun.taylorTesting(func,2);

% We expect all elements of ORDER1 to be close to 1. Since the methods being
% tested in this case are all linear, ORDER2 will be noise. However, since
% the methods are indeed linear, we should expect nDiff2 to have values all
% close to machine epsilon, which we can use to check for the correctness of
% the derivative computed.
pass(2) = ( (max(abs(order1 - 1)) < tolOrder) && ...
    (max(abs(nDiff2)) < tolDiff) );


% Check that we received the correct linearity information
pass(3) = ( lin == 1 );

end
