% Test file for ADCHEBFUN times

function pass = test_times

% Functions to be tested
funcList = {@times};

% Initialise pass vector
pass = zeros(2*length(funcList), 5);

% Steps to be taken in Taylor testing
numSteps = 6;

% Tolerance for Taylor testing
tolOrder = 1e-2;
tolDiff = 1e-14;

% Loop through functions
for funcCounter = 1:length(funcList)
    % Temporarily store the function handle we're working with:
    func = funcList{funcCounter};
    
    % Compare values
    err = adchebfun.valueTestingBinary(func);
    
    % Confirm that the error returned is zero
    pass(2*(funcCounter-1) + 1, :) = ( err == 0 );
    
    % Taylor testing
    [order1, order2, nDiff2] = adchebfun.taylorTestingBinary(func, numSteps);
    
    % We expect all elements of ORDER1 to be close to 1. Since the methods being
    % tested in this case are all linear, ORDER2 will be noise. However, since
    % the methods are indeed linear, we should expect nDiff2 to have values all
    % close to machine epsilon, which we can use to check for the correctness of
    % the derivative computed.
    
    % Check whether we get the expected results for linear and nonlinear
    % operations.
    linearOpResults = ( ((max(abs(order1(:, 2:end))) - 1) < tolOrder) & ...
        (max(abs(nDiff2(:, 2:end))) < tolDiff) );
    
    nonlinearOpResults = ( (max(abs(order1(:, 1) - 1)) < tolOrder) && ...
        (max(abs(order2(:, 1) - 2)) < tolOrder) );
    
    % Concatenate results
    pass(2*(funcCounter-1) + 2, :) = [nonlinearOpResults, linearOpResults];
    
end

end