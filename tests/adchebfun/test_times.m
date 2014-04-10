% Test file for ADCHEBFUN times

function pass = test_times

% Initialise pass vector
pass = zeros(3, 5);

% Tolerance for Taylor testing
tolOrder = 1e-2;
tolDiff = 1e-14;

% The function handle we're working with:
func = @times;

% Compare values
[err, lin] = adchebfun.valueTestingBinary(func);

% Confirm that the error returned is zero
pass(1, :) = ( err == 0 );

% Taylor testing
[order1, order2, nDiff2] = adchebfun.taylorTestingBinary(func);

% We expect all elements of ORDER1 to be close to 1. Since depending on the
% combinations of variables (ADCHEBFUN, CHEBFUN and SCALAR), we should expect
% some of ORDER2 to be noise (for linear operations), in which case NDIFF2
% should have values close to machine epsilon. For the nonlinear operations,
% ORDER2 will be close having all elements value close to 2.

% Check whether we get the expected results for linear and nonlinear
% operations.
linearOpResults = ( ((max(abs(order1(:, 2:end))) - 1) < tolOrder) & ...
    (max(abs(nDiff2(:, 2:end))) < tolDiff) );

nonlinearOpResults = ( (max(abs(order1(:, 1) - 1)) < tolOrder) && ...
    (max(abs(order2(:, 1) - 2)) < tolOrder) );

% Concatenate results
pass(2, :) = [nonlinearOpResults, linearOpResults];

% Linearity checking
pass(3, :) = ( lin == [0, 1, 1, 1, 1] );
end