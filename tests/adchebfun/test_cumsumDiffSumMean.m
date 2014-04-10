% Test file for ADCHEBFUN CUMSUM, DIFF, MEAN and SUM

function pass = test_cumsumDiffSumMean

% List of trigonometric functions to test.
funcList = {@diff, @(u)diff(u, 2), @(u)diff(u, 4), ...
                 @sum, @(u) sum(u, -.25, .8), ...
                 @cumsum, @(u)cumsum(u,2), ...
                 @mean};

% Tolerance for Taylor testing
tolOrder = 1e-2;
tolDiff = 1e-12;
% Initialise vector with pass information
pass = zeros(2, numel(funcList));

% Do the tests.
for k = 1:numel(funcList)
    % Call the valueTesting method, which also returns linearity information
    [err, lin] = adchebfun.valueTesting(funcList{k});
    
    % First, check that the computed function values match what we expect
    pass(1, k) = ( err == 0 );
    
    % Call the taylorTesting method
    [order1, order2, nDiff2] = adchebfun.taylorTesting(funcList{k});
    
    % We expect all elements of ORDER1 to be close to 1. Since the methods being
    % tested in this case are all linear, ORDER2 will be noise. However, since
    % the methods are indeed linear, we should expect nDiff2 to have values all
    % close to machine epsilon, which we can use to check for the correctness of
    % the derivative computed.
    pass(2,k) = ( (max(abs(order1 - 1)) < tolOrder) && ...
        (max(abs(nDiff2)) < tolDiff) );
    
    % Check that we received the correct linearity information
    pass(3, k) = ( lin == 1 );
end

end