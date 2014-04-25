% Test file for ADCHEBFUN trigonometric and related functions, part 3/4

function pass = test_trig3

% List of trigonometric functions to test.
funcList = {@atand, @atanh, @cos, @cosd, @cosh, @cot, @cotd, @coth, ...
    @csc, @cscd};

% Call the ADCHEBFUN testUnary() method to do the tests.
pass = adchebfun.testUnary(funcList);

end