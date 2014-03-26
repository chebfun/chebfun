% Test file for ADCHEBFUN times

function pass = test_power

% Initialise pass vector
pass = zeros(3, 5);

% Tolerance for Taylor testing
tolOrder = 1e-2;
tolDiff = 1e-14;

% The function handle we're working with:
func = @power;

% Compare values
[err, lin] = adchebfun.valueTestingBinary(func);

% Confirm that the error returned is zero
pass(1, :) = ( err < tolDiff );

% Taylor testing
[order1, order2] = adchebfun.taylorTestingBinary(func);

% We expect all elements of ORDER1 to be close to 1. Since the methods being
% tested in this all nonlinear, ORDER2 should have entries close to 2.

% Check whether we get the expected results for linear and nonlinear
% operations.
pass(2,:) = ( (max(abs(order1 - 1)) < tolOrder) & ...
              (max(abs(order2 - 2)) < tolOrder) );

% Linearity checking. All operations should be detected as nonlinear.
pass(3,:) = ( lin == 0 );
end