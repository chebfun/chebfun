% Test file for ADCHEBFUN times

function pass = test_times

% Initialise pass vector
pass = zeros(2, 5);

% Tolerance for Taylor testing
tolOrder = 1e-2;
tolDiff = 1e-14;

% The function handle we're working with:
func = @times;

% Compare values
err = adchebfun.valueTestingBinary(func);

% Confirm that the error returned is zero
pass(1, :) = ( err == 0 );

% Taylor testing
[order1, order2, nDiff2] = adchebfun.taylorTestingBinary(func);

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
pass(2, :) = [nonlinearOpResults, linearOpResults];


end