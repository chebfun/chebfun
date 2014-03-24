% Test file for ADCHEBFUN times

function pass = test_rdivide

% Initialise pass vector
pass = zeros(2, 5);

% Tolerance for Taylor testing
tolOrder = 1e-2;
tolDiff = 1e-14;

% The function handle we're working with:
func = @rdivide;

% Compare values
[err, lin] = adchebfun.valueTestingBinary(func);

% Confirm that the error returned is zero
pass(1, :) = ( err < tolDiff );

% Taylor testing
[order1, order2, nDiff2] = adchebfun.taylorTestingBinary(func);

% We expect all elements of ORDER1 to be close to 1. Since depending on the
% combinations of variables (ADCHEBFUN, CHEBFUN and SCALAR), we should expect
% some of ORDER2 to be noise (for linear operations), in which case NDIFF2
% should have values close to machine epsilon. For the nonlinear operations,
% ORDER2 will be close having all elements value close to 2.

% Check whether we get the expected results for linear and nonlinear
% operations.
linearOpResults = ( ((max(abs(order1)) - 1) < tolOrder) & ...
    (max(abs(nDiff2)) < tolDiff) );

nonlinearOpResults = ( (max(abs(order1 - 1)) < tolOrder) & ...
    (max(abs(order2 - 2)) < tolOrder) );

% Cherry pick the results from the vectors above depending on whether the
% operations were linear or nonlinear.
pass(2,1) = nonlinearOpResults(1);      % ADCHEBFUN and ADCHEBFUN
pass(2,2) = linearOpResults(2);         % ADCHEBFUN and CHEBFUN
pass(2,3) = nonlinearOpResults(3);      % CHEBFUN and ADCHEBFUN
pass(2,4) = linearOpResults(4);         % ADCHEBFUN and SCALAR
pass(2,5) = nonlinearOpResults(5);      % SCALAR and ADCHEBFUN

% Linearity checking
pass(3, :) = ( lin == [0, 1, 0, 1, 0] );

end