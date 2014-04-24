% Test file for ADCHEBFUN trigonometric and related functions, part 2/4

function pass = test_trig2

% List of trigonometric functions to test.
funcList = {@acsch, @asec, @asecd, @asech, @asin, @asind, @asinh, @atan};

% Call the ADCHEBFUN testUnary() method to do the tests.
pass = adchebfun.testUnary(funcList);

end