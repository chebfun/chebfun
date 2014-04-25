% Test file for ADCHEBFUN CUMPROD() and PROD()

function pass = test_cumprodProd

% List of trigonometric functions to test.
funcList = {@cumprod, @prod};

% Call the ADCHEBFUN testUnary() method to do the tests.
pass = adchebfun.testUnary(funcList);

end