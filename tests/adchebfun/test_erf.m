% Test file for ADCHEBFUN error and related functions

function pass = test_erf

% List of functions related to the error function to test.
funcList = {@erf, @erfc, @erfcinv, @erfcx, @erfinv};

% Call the ADCHEBFUN testUnary() method to do the tests.
pass = adchebfun.testUnary(funcList);

end