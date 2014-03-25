% Test file for ADCHEBFUN exp, log related functions

function pass = test_expLog

% List of trigonometric functions to test.
funcList = {@exp, @expm1, @log, @log10, @log1p, @log2};

% Call the ADCHEBFUN testUnary() method to do the tests.
pass = adchebfun.testUnary(funcList);

end