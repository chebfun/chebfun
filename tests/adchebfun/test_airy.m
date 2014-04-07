% Test file for ADCHEBFUN Airy function

function pass = test_airy

% List of airy functions to test, with different number of inputs
airy2 = @(f) airy(2, f);
funcList = {@airy, airy2};

% Call the ADCHEBFUN testUnary() method to do the tests.
pass = adchebfun.testUnary(funcList);

end
