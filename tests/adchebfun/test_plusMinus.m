% Test file for ADCHEBFUN plus and minus

function pass = test_plusMinus

% Functions to be tested
funcList = {@plus, @minus};

% Initialise pass vector
pass = zeros(3*numel(funcList), 5);

% Tolerance for Taylor testing
tolOrder = 1e-2;
tolDiff = 1e-14;

% Loop through functions
for funcCounter = 1:length(funcList)
    % Temporarily store the function handle we're working with:
    func = funcList{funcCounter};
    
    % Compare values
    [err, lin] = adchebfun.valueTestingBinary(func);
    
    % Confirm that the error returned is zero
    pass(3*(funcCounter-1) + 1, :) = ( err == 0 );
    
    % Taylor testing
    [order1, order2, nDiff2] = adchebfun.taylorTestingBinary(func);
    
    % We expect all elements of ORDER1 to be close to 1. Since the methods being
    % tested in this case are all linear, ORDER2 will be noise. However, since
    % the methods are indeed linear, we should expect nDiff2 to have values all
    % close to machine epsilon, which we can use to check for the correctness of
    % the derivative computed.
    pass(3*(funcCounter-1) + 2, :) = ( (max(abs(order1 - 1)) < tolOrder) & ...
        (max(abs(nDiff2)) < tolDiff) );
    
    % Linearity checking. All operations should be detected as linear.
    pass(3*(funcCounter-1) + 3, :) = ( lin == 1);
    
end

end