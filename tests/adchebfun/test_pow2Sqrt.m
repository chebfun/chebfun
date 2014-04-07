% Test file for ADCHEBFUN POW2() and SQRT()

function pass = test_pow2Sqrt

% List of trigonometric functions to test.
funcList = {@pow2, @sqrt};

% Call the ADCHEBFUN testUnary() method to do the tests.
pass = adchebfun.testUnary(funcList);

end