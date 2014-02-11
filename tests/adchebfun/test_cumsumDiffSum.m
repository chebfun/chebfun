% Test file for ADCHEBFUN CUMSUM, DIFF, and SUM

function pass = test_cumsumDiffSum

% List of trigonometric functions to test.
diffFunctions = {@diff, @(u)diff(u,2), @(u)diff(u,4), ...
                 @sum, ...
                 @cumsum, @(u)cumsum(u,2)};

% How many iterations we want to in the Taylor testing
numSteps = 4;

% Tolerance for Taylor testing
tolOrder = 1e-2;
tolDiff = 1e-12;
% Initialise vector with pass information
pass = zeros(2, numel(diffFunctions));

% Do the tests.
for k = 1:numel(diffFunctions)
    % First, check that the computed function values match what we expect
    pass(1,k) = ( adchebfun.valueTesting(diffFunctions{k}) == 0 );
    
    % Call the taylorTesting method
    [order1, order2, nDiff2] = adchebfun.taylorTesting(diffFunctions{k}, numSteps);
    
    % We expect all elements of ORDER1 to be close to 1. Since the methods being
    % tested in this case are all linear, ORDER2 will be noise. However, since
    % the methods are indeed linear, we should expect nDiff2 to have values all
    % close to machine epsilon, which we can use to check for the correctness of
    % the derivative computed.
    pass(2,k) = ( max(abs(order1 - 1)) < tolOrder && ...
        max(abs(nDiff2)) < tolDiff);
end

end