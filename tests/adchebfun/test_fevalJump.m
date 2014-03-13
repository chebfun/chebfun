% Test file for ADCHEBFUN feval and jump functions

function pass = test_fevalJump

% List of airy functions to test, with different number of inputs

fevalFun = @(u) feval(u, .3);
jumpFun = @(u) jump(u, .5, -1);
% List of trigonometric functions to test.
funcList = {fevalFun, jumpFun};

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
    %
    % However, the jump method will run into issues of dividing 0 by 0. This
    % will result in NaNs, which can only occur if it gets the derivative
    % correct. So check instead that nDiff2 has all zero entries.
    if ( k ~= 2)
        pass(2,k) = ( (max(abs(order1 - 1)) < tolOrder) && ...
        (max(abs(nDiff2)) < tolDiff) );
    else
        pass(2,2) = ~any(nDiff2);
    end
    
    % Check that we received the correct linearity information
    pass(3, k) = ( lin == 1 );
end

end