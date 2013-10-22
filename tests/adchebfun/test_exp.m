% Test file for adchebfun/sin

function pass = test_exp

% First, check that the computed function values match what we expect
pass(1) = ( valueTesting(@exp) == 0 );

% Call the taylorTesting method
[order1, order2] = taylorTesting(@exp,7);

% We expect all elements of ORDER1 to be close to 1, and of ORDER2 to be close
% to 2.
pass(2) = ( max(abs(order1-1)) < .01 ) & ( max(abs(order2-2)) < .01 );

end

